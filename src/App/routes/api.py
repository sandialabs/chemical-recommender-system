# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pubchempy as pcp
from flask import Blueprint, jsonify, Response, current_app, stream_with_context, request
import queue
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem.Draw import rdMolDraw2D
import io
import base64
import json
import time
import logging
import matplotlib
matplotlib.use('Agg')  # Force headless backend for server
from Comparison.utils.gen import parseQuery as pq
from utils.progress_logger import ProgressLogger

# Define Blueprint and logger at the very top, before any route decorators
api_bp = Blueprint("api", __name__)
logger = logging.getLogger(__name__)

def get_progress_queue(job_id):
    """Get the progress queue from the ProgressLogger system."""
    job_id_str = str(job_id)
    
    # Use the ProgressLogger's queue storage system instead of a separate one
    with ProgressLogger._lock:
        if job_id_str in ProgressLogger._progress_queues:
            logger.debug(f"Found existing queue for job {job_id_str}")
            return ProgressLogger._progress_queues[job_id_str]
        else:
            # Create a new queue and store it in ProgressLogger's system
            import queue
            new_queue = queue.Queue()
            ProgressLogger._progress_queues[job_id_str] = new_queue
            logger.debug(f"Created new queue for job {job_id_str}")
            return new_queue

@api_bp.route("/api/similarity_map")
def similarity_map():
    """
    Returns a similarity map image (PNG) comparing the given CID to the query molecule (reference).
    Query molecule is taken from current_app.state.queryval.
    Target molecule is given by ?cid=...
    """
    state = current_app.state
    queryval = getattr(state, "queryval", None)
    target_cid = request.args.get("cid", None)
    if not queryval or not target_cid:
        return jsonify({"error": "Missing query or cid"}), 400

    # Get SMILES for query and target
    import traceback
    try:
        # Query molecule
        from Comparison.utils.gen import get_compound_properties
        if queryval.isdigit():
            query_smiles, _ = get_compound_properties(int(queryval))
        else:
            # Try to resolve as name or SMILES
            try:
                query_cid = pcp.get_cids(queryval, "name", list_return="flat")[0]
                query_smiles, _ = get_compound_properties(query_cid)
            except Exception as e:
                logger.error(f"Could not resolve queryval '{queryval}' as name: {e}")
                query_smiles = queryval

        # Target molecule
        try:
            target_smiles, _ = get_compound_properties(int(target_cid))
        except Exception as e:
            logger.error(f"Could not resolve target_cid '{target_cid}': {e}")
            return jsonify({"error": f"Could not resolve target_cid '{target_cid}': {e}"}), 400

        logger.info(f"Query SMILES: {query_smiles}, Target SMILES: {target_smiles}")
        refmol = Chem.MolFromSmiles(query_smiles)
        mol = Chem.MolFromSmiles(target_smiles)
        if not refmol or not mol:
            logger.error(f"Could not parse molecules: query_smiles={query_smiles}, target_smiles={target_smiles}")
            return jsonify({"error": f"Could not parse molecules: query_smiles={query_smiles}, target_smiles={target_smiles}"}), 400

        # Create a 2D drawing object
        d2d = Draw.MolDraw2DCairo(300, 300)

        # Generate the similarity map
        try:
            _, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(
                refmol, mol,
                fpFunction = lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv'),
                draw2d=d2d,
            )
            d2d.FinishDrawing()
        except Exception as e:
            logger.error(f"Error in GetSimilarityMapForFingerprint: {e}")
            # fallback: draw plain molecule
            d2d.DrawMolecule(mol)
            d2d.FinishDrawing()
        img_bytes = d2d.GetDrawingText()

        # Return as base64 data URL for easy frontend use
        img_b64 = base64.b64encode(img_bytes).decode('utf-8')
        return jsonify({"image": "data:image/png;base64," + img_b64})
    except Exception as e:
        tb = traceback.format_exc()
        logger.error(f"Error generating similarity map: {e}\n{tb}")
        return jsonify({"error": str(e), "traceback": tb}), 500
    
@api_bp.route("/progress_stream")
def progress_stream():
    job_id = str(request.args.get("job_id", "default"))
    logger.info(f"Progress stream requested for job_id: {job_id}")
    
    # Make sure the queue exists
    q = get_progress_queue(job_id)

    def event_stream():
        # Set connection timeout longer
        timeout = 30
        last_ping = time.time()
        ping_interval = 5  # Send a ping every 5 seconds
        
        try:
            # Start streaming progress messages
            while True:
                try:
                    # Check if we need to send a keep-alive ping
                    current_time = time.time()
                    if current_time - last_ping > ping_interval:
                        yield ": ping\n\n"
                        last_ping = current_time
                    
                    # Try to get a message with a shorter timeout to allow for pings
                    try:
                        message = q.get(block=True, timeout=1)
                        logger.debug(f"Sending progress update for job {job_id}: {message}")
                        yield f"data: {message}\n\n"
                    except queue.Empty:
                        # No message yet, just continue the loop for pings
                        continue
                    
                except Exception as e:
                    # Log any errors but try to continue
                    logger.error(f"Error in event stream loop for {job_id}: {str(e)}")
                    yield f"data: {json.dumps({'status': 'Warning', 'detail': f'Stream error: {str(e)}'})}\n\n"
                    time.sleep(1)  # Prevent tight loop if there's an error
        
        except GeneratorExit:
            # Client disconnected
            logger.info(f"Client disconnected from SSE stream for job {job_id}")
        except Exception as e:
            logger.error(f"Fatal error in event stream for job {job_id}: {str(e)}")
            yield f"data: {json.dumps({'status': 'Error', 'detail': f'Stream terminated: {str(e)}'})}\n\n"

    # Set response headers for SSE
    response = Response(
        stream_with_context(event_stream()),
        mimetype="text/event-stream"
    )
    # Add additional headers to prevent buffering
    response.headers['Cache-Control'] = 'no-cache, no-transform'
    response.headers['X-Accel-Buffering'] = 'no'  # For Nginx
    return response

@api_bp.route("/api/query_image")
def query_image():
    print("started")
    state = current_app.state
    query_smiles = state.queryval

    mol = Chem.MolFromSmiles(query_smiles)

    img_width = 300
    img_height = 300
    drawer = rdMolDraw2D.MolDraw2DCairo(img_width, img_height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    img_bytes = drawer.GetDrawingText()

    # Encode image bytes as base64 for JSON transport
    img_b64 = base64.b64encode(img_bytes).decode('utf-8')
    print("returning")
    return jsonify({"image": "data:image/png;base64," + img_b64})