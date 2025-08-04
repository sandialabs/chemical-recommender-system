# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

from flask import Blueprint, render_template, request, jsonify, current_app
from Comparison.Controller import BatchRun as br
import threading
import uuid
import os
import logging

BATCH_JOBS_DIR = "src/Comparison/LocalIO/jobs"

batch_bp = Blueprint("batch", __name__)
logger = logging.getLogger(__name__)


def batch_processing(app, batch_text, job_id, job_dir):
    """
    Function to handle batch processing in a separate process.
    Now writes all outputs to a per-job directory and exposes them via static symlink.

    Args:
        app (Flask): The Flask application context.
        batch_text (str): The input text for batch processing.
        batch_status (multiprocessing.Value): A shared value to indicate the status of the batch processing.
            - "loading": Indicates that the batch processing is currently running.
            - "true": Indicates that the batch processing completed successfully.
            - "false": Indicates that there was an error during batch processing.
        job_id (str): The unique job ID for this batch job.
        job_dir (str): The directory where all outputs for this job will be written.
    """
    with app.app_context():
        # Import here to avoid circular imports
        from App.routes.api import get_progress_queue
        
        try:
            # Get the progress queue for this job
            progress_queue = get_progress_queue(job_id)
            
            # Send initial progress message
            import json
            try:
                initial_msg = json.dumps({
                    "status": "Initializing", 
                    "detail": "Setting up batch processing environment..."
                })
                if progress_queue:
                    progress_queue.put_nowait(initial_msg)
            except Exception as e:
                logger.debug(f"Could not send initial progress message: {e}")
            
            logger.info(f"Starting BatchRun with output path: {os.path.join(job_dir, 'batch')}")
            
            # Pass job_dir to BatchRun so all outputs are written there
            br(batch_text, output=os.path.join(job_dir, "batch"), job_id=job_id, progress_queue=progress_queue)
            
            # BatchRun will send the final completion message, so we don't need to send another here
            logger.info("Batch processing thread completed - BatchRun handles final status message")
            logger.info(f"Checking for output files in {job_dir}")
            
            # List the files created in the job directory
            if os.path.exists(job_dir):
                files = os.listdir(job_dir)
                logger.info(f"Files created in job directory: {files}")
            else:
                logger.error(f"Job directory does not exist: {job_dir}")
            
            # After job completes, symlink/copy job_dir to static/jobs/<job_id>
            static_jobs_dir = os.path.join(app.root_path, "static", "jobs")
            os.makedirs(static_jobs_dir, exist_ok=True)
            static_job_dir = os.path.join(static_jobs_dir, job_id)
            logger.info(f"Creating static symlink from {job_dir} to {static_job_dir}")
            
            if not os.path.exists(static_job_dir):
                try:
                    os.symlink(os.path.abspath(job_dir), static_job_dir)
                    logger.info(f"Symlink created successfully")
                except Exception as e:
                    logger.warning(f"Symlink failed, trying copy: {e}")
                    # If symlink fails (e.g., on Windows), fallback to copytree
                    import shutil
                    shutil.copytree(job_dir, static_job_dir)
                    logger.info(f"Copy completed successfully")
            
            # Update the shared batch status through the app state
            # Note: we need to be very careful about updating the shared state from a thread
            logger.info(f"Attempting to update batch_status to 'true' for job {job_id}")
            try:
                with app.app_context():
                    from flask import current_app
                    current_app.state.batch_status = "true"
                    logger.info(f"Successfully updated batch_status to: {current_app.state.batch_status}")
            except Exception as state_error:
                logger.error(f"Failed to update state: {state_error}")
                # Fallback: try to update the original app state reference
                app.state.batch_status = "true"
                logger.info(f"Fallback: set app.state.batch_status to: {app.state.batch_status}")
            
            logger.info(f"Batch processing completed successfully for job {job_id}, status set to: true")
            
        except Exception as e:
            # Update the shared batch status on error through the app state
            logger.error(f"Batch processing failed for job {job_id}: {e}")
            try:
                with app.app_context():
                    from flask import current_app
                    current_app.state.batch_status = "false"
            except Exception as state_error:
                logger.error(f"Failed to update error state: {state_error}")
                app.state.batch_status = "false"
            
            # Send error message
            try:
                error_msg = json.dumps({
                    "status": "Error", 
                    "detail": f"Batch processing failed: {str(e)}"
                })
                if progress_queue:
                    progress_queue.put_nowait(error_msg)
            except Exception as progress_error:
                logger.debug(f"Could not send error progress message: {progress_error}")
            
            logger.error(f"Error during batch processing: {e}")


@batch_bp.route("/batch", methods=["GET", "POST"])
def batch():
    """
    Route to handle batch processing requests.

    Methods:
        GET: Renders the batch processing page.
        POST: Starts a new batch processing task with the provided input text.

    Returns:
        Rendered template for the batch processing page with the current status.

    Args:
        None

    Request Form Data (POST):
        batchinput (str): The input text for batch processing.

    Context Variables:
        state (object): The current application state, which includes:
            - batch_status (multiprocessing.Value): The current status of the batch processing.
            - batch_process (multiprocessing.Process): The current batch processing process.
    """
    state = current_app.state
    state.reset_vals()
    if request.method == "POST":
        batch_text = request.form.get("batchinput")
        # Generate a unique job ID for this batch job
        job_id = str(uuid.uuid4())
        job_dir = os.path.join(BATCH_JOBS_DIR, job_id)
        os.makedirs(job_dir, exist_ok=True)
        state.batch_status = "loading"
        state.batch_job_id = job_id  # Track the current job ID
        if state.batch_process and state.batch_process.is_alive():
            state.batch_process.join(timeout=1)
        # Pass job_id and job_dir to the batch thread
        state.batch_process = threading.Thread(
            target=batch_processing,
            args=(current_app._get_current_object(), batch_text, job_id, job_dir),
            daemon=True
        )
        state.batch_process.start()
        # Pass job_id to the frontend for polling
        return render_template("batch.html", ready=state.batch_status, job_id=job_id)
    return render_template("batch.html", ready=state.batch_status)


@batch_bp.route("/batch_status", methods=["GET"])
def batch_status():
    """
    Route to get the current status of the batch processing.

    Methods:
        GET: Returns the current status of the batch processing.

    Returns:
        JSON response with the current status of the batch processing.

    Args:
        None

    Context Variables:
        state (object): The current application state, which includes:
            - batch_status (multiprocessing.Value): The current status of the batch processing.
    """
    state = current_app.state
    logger.info(f"Batch status request - status: {state.batch_status}, job_id: {state.batch_job_id}")
    return jsonify({
        "ready": state.batch_status, 
        "job_id": state.batch_job_id
    })


@batch_bp.route("/batch_restart", methods=["POST"])
def batch_restart():
    """
    Route to restart the batch processing.

    Methods:
        POST: Terminates any ongoing batch processing and resets the state.

    Returns:
        JSON response with the current status of the batch processing.

    Args:
        None

    Context Variables:
        state (object): The current application state, which includes:
            - batch_status (multiprocessing.Value): The current status of the batch processing.
            - batch_process (multiprocessing.Process): The current batch processing process.
    """
    state = current_app.state
    if state.batch_process and state.batch_process.is_alive():
        state.batch_process.join(timeout=1)
    state.batch_process = None
    state.reset_vals()
    return jsonify({"ready": state.batch_status})
