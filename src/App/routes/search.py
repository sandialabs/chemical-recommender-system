# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

from flask import (
    Blueprint,
    render_template,
    request,
    current_app,
    jsonify,
)
import uuid
import json
import threading
import time
import logging

search_bp = Blueprint("search", __name__)
logger = logging.getLogger(__name__)

def search_processing(
    app,
    query,
    smarts,
    smarts_num,
    finnum,
    tarray,
    incEle,    
    disallow_isotopes,
    result_queue,
    job_id="default",
):
    """
    Runs the search processing function.

    Args:
    - app: The Flask application instance.
    - query (str): The query for search processing.
    - smarts (str): SMARTS substructure search.
    - smarts_num (str): Number of SMARTS matches.
    - finnum (int): Number of final candidates.
    - tarray (list): Thermophysical properties array.
    - incEle (list): Included elements.
    - search_status: A shared value to update the search status.
    - captured_output: A shared StringIO object to capture the standard output.
    - result_queue: A multiprocessing.Queue to pass results back to the main thread.
    - job_id (str): The job ID for the search processing.

    This function updates the search status to "true" after processing.
    """
    with app.app_context():
        # Import here to avoid circular imports
        from Comparison.Controller import SingleRun
        from App.routes.api import get_progress_queue
        
        logger = logging.getLogger(__name__)
        
        try:
            # Get the progress queue for this job
            progress_queue = get_progress_queue(job_id)
            
            # Send initial progress message
            logger.debug(f"Starting search processing for job {job_id}, query: {query}")
            progress_queue.put(json.dumps({
                "status": "Starting", 
                "detail": f"Initializing search for query: {query}"
            }))
            
            # Use SingleRun which has proper progress reporting
            weights = [1, 1, 1, 1, 1]  # Default weights
            containers = []  # Default containers
            
            logger.debug(f"About to call SingleRun for job {job_id}")
            
            # Run SingleRun which handles all the processing, file generation, AND returns results
            results_data = SingleRun(
                query, finnum, tarray, incEle, smarts, smarts_num, 
                weights, containers, job_id, progress_queue, disallow_isotopes
            )
            
            logger.debug(f"SingleRun completed for job {job_id}")
            
            # Check if we got valid results
            if results_data is None:
                logger.warning(f"SingleRun returned None for job {job_id}")
                # Error occurred, SingleRun already sent error message to progress queue
                result_queue.put({
                    "status": "false", 
                    "error": "Search failed - see progress messages for details", 
                    "job_id": job_id
                })
                return
            
            # Unpack the results - now includes opera_failed
            if len(results_data) == 4:
                results, query_models, subfailed, opera_failed = results_data

            else:
                logger.error(f"Unexpected results_data length: {len(results_data)}")
                result_queue.put({
                    "status": "false", 
                    "error": "Unexpected results format", 
                    "job_id": job_id
                })
                return
            
            logger.debug(f"SingleRun returned results with length: {len(results) if results else 0}")
            logger.debug(f"Results type: {type(results)}")
            logger.debug(f"OPERA failed: {opera_failed}")
            if results and len(results) > 0:
                logger.debug(f"First result in search processing: {results[0] if len(results) > 0 else 'None'}")
            
            # Add OPERA warning to progress if needed
            if opera_failed:
                progress_queue.put(json.dumps({
                    "status": "Warning", 
                    "detail": "Warning: Property predictions failed - results have limited accuracy for thermal/toxicity comparisons"
                }))
            
            # Send completion message with redirect URL
            progress_queue.put(json.dumps({
                "status": "Finished", 
                "detail": "Search completed successfully. Loading your results..." + (" (with property calculation warnings)" if opera_failed else ""),
                "redirect_url": f"/results?job_id={job_id}"
            }))
            logger.debug(f"Sent 'Finished' status to progress queue for job {job_id}")
            
            # Put results in the result queue for the main thread
            result_data = {
                "queryval": query,
                "results": results,
                "query_models": query_models,
                "subfailed": subfailed,
                "opera_failed": opera_failed, 
                "status": "true",
                "from_ui": "true",
                "job_id": job_id,
            }
            logger.debug(f"About to put in result_queue - queryval: '{query}'")
            logger.debug(f"query type: {type(query)}, query length: {len(query) if query else 0}")
            result_queue.put(result_data)
            
            logger.debug(f"Search processing completed successfully for job {job_id}")
            logger.debug(f"Results put in result_queue: status=true, results_length={len(results) if results else 0}")
            logger.debug(f"Result queue size after put: {result_queue.qsize()}")
                
        except Exception as e:
            # Send error message to progress queue
            try:
                progress_queue = get_progress_queue(job_id)
                progress_queue.put(json.dumps({
                    "status": "Error", 
                    "detail": f"Search failed: {str(e)}"
                }))
            except:
                pass  # Don't let progress messaging errors mask the original error
            
            result_queue.put({
                "status": "false", 
                "error": str(e), 
                "job_id": job_id
            })
            logger.error(f"Error during search processing for job {job_id}: {e}")
            import traceback
            traceback.print_exc()


@search_bp.route("/search", methods=["GET", "POST"])
def search():
    """
    Handles the search functionality for the application.

    This function supports both GET and POST requests:
    - GET: Renders the search page.
    - POST: Processes the search query and related parameters, performs the comparison, and updates the status.

    Inputs:
    - None directly (uses request.form to get input data from the POST request)

    Application State:
    - state.queryval: The query value provided by the user.
    - state.finnum: The number of candidate results to be returned by the run.
    - state.tarray: An array of boolean values indicating which thermophysical properties are to be calculated.
    - state.weights: An array of weights corresponding to the importance of each property in tarray.
    - state.params: Parameters for the comparison function.
    - state.results: The results of the comparison function.
    - state.query_models: Models used in the query.
    - state.subfailed: A flag indicating if substructure searching failed.

    Returns:
    - GET: Renders the 'search.html' template.
    - POST: Renders the 'search.html' template with the status.
    """
    logger = logging.getLogger(__name__)
    state = current_app.state
    state.reset_vals()
    
    if request.method == "POST":
        # Get job_id from form data (sent by frontend) or generate new one
        job_id = request.form.get("job_id")
        if not job_id:
            job_id = f"search_{uuid.uuid4().hex}"
        state.single_job_id = job_id
        
        # Retrieve form data
        query = request.form.get("query")
        state.queryval = query
        smarts = request.form.get("smarts")
        if smarts == "":
            smarts = None
        smarts_num = request.form.get("smarts_num")
        if smarts_num == "":
            smarts_num = None
        state.finnum = int(request.form.get("finals"))

        # Retrieve and process thermophysical property selections
        MP = request.form.get("MP") == "on"
        BP = request.form.get("BP") == "on"
        logP = request.form.get("logP") == "on"
        Hlaw = request.form.get("HLaw") == "on"
        VP = request.form.get("VP") == "on"

        incEle = request.form.get("IncEle") == "on"
        incEle2 = request.form.get("IncEle2")
        disallow_isotopes = request.form.get("disallow_isotopes") == "on"

        if not incEle:
            incEle = [element.strip() for element in incEle2.split(",")]

        # Update application state with the selected properties and default weights
        state.tarray = [MP, BP, logP, Hlaw, VP]
        state.weights = [1, 1, 1, 1, 1]
        state.params = [
            query,
            state.finnum,
            state.tarray,
            incEle,
            smarts,
            smarts_num,
            disallow_isotopes,
        ]


        if query:
            state.search_status = "loading"
            logger.debug(f"Starting search processing for job {job_id} at {time.time()}")
            if state.search_process and state.search_process.is_alive():
                state.search_process.join(timeout=1)
            # Create a queue to receive results from the worker thread
            import queue
            result_queue = queue.Queue()
            # Start the search processing in a separate thread
            state.search_process = threading.Thread(
                target=search_processing,
                args=(
                    current_app._get_current_object(),
                    query,
                    smarts,
                    smarts_num,
                    state.finnum,
                    state.tarray,
                    incEle,
                    disallow_isotopes,
                    result_queue,
                    job_id,  # Pass job_id to search_processing
                ),
                daemon=True
            )
            state.search_process.start()
            logger.debug(f"Background thread started for job {job_id} at {time.time()}")
            # Store the result queue in the state for later access
            state.result_queue = result_queue
            
            if request.headers.get('X-Requested-With') == 'XMLHttpRequest' or request.headers.get('Content-Type') == 'application/json':
                return jsonify({
                    "status": "success",
                    "job_id": job_id,
                    "ready": state.search_status
                })
            else:
                # Return the search template with loading status and job_id for regular form submission
                logger.debug(f"Returning search template for job {job_id} at {time.time()}")
                return render_template("search.html", ready=state.search_status, job_id=job_id)
    else:
        # For GET requests, generate a default job_id
        job_id = request.args.get("job_id", f"search_{uuid.uuid4().hex}")
        state.single_job_id = job_id

    # Render the search.html template for GET requests with job_id
    return render_template("search.html", ready=state.search_status, job_id=job_id)


@search_bp.route("/search_status", methods=["GET"])
def search_status():
    """
    Returns the current status of the search processing.

    This endpoint is used by the client to poll the server for the status of the search processing.

    Returns:
    - JSON object with the 'ready' status.
    """
    state = current_app.state

    # Check if there are results in the queue and update the state
    if state is not None and hasattr(state, "result_queue") and state.result_queue is not None and not state.result_queue.empty():
        logger.debug("Found results in result_queue, processing...")
        result = state.result_queue.get()
        logger.debug(f"Retrieved result from queue: status={result.get('status')}, job_id={result.get('job_id')}")
        
        if result["status"] == "true":
            logger.debug("Setting state values from successful result")
            state.queryval = result["queryval"]
            state.results = result["results"]
            state.query_models = result["query_models"]
            state.subfailed = result["subfailed"]
            state.opera_failed = result.get("opera_failed", False)  # Add OPERA failure status to state
            state.search_status = "true"
            
            logger.debug(f"Updated state - queryval: '{state.queryval}', results count: {len(state.results) if state.results else 0}, opera_failed: {state.opera_failed}")
        else:
            logger.warning(f"Result status was not 'true', resetting state. Error: {result.get('error')}")
            state.reset_vals()
            logger.error(f"Error during search processing: {result.get('error')}")
        state.search_process = None

    return jsonify({"ready": state.search_status if state is not None else "false"})


