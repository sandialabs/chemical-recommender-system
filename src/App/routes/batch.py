# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

from flask import Blueprint, render_template, request, jsonify, current_app
from Comparison.Controller import BatchRun as br
import multiprocessing

batch_bp = Blueprint("batch", __name__)


def batch_processing(app, batch_text, batch_status):
    """
    Function to handle batch processing in a separate process.

    Args:
        app (Flask): The Flask application context.
        batch_text (str): The input text for batch processing.
        batch_status (multiprocessing.Value): A shared value to indicate the status of the batch processing.
            - "loading": Indicates that the batch processing is currently running.
            - "true": Indicates that the batch processing completed successfully.
            - "false": Indicates that there was an error during batch processing.
    """
    with app.app_context():
        try:
            br(batch_text)
            batch_status.value = "true"
        except Exception as e:
            batch_status.value = "false"
            print(f"Error during batch processing: {e}")


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
        state.batch_status.value = "loading"
        if state.batch_process and state.batch_process.is_alive():
            state.batch_process.terminate()
            state.batch_process.join()
        state.batch_process = multiprocessing.Process(
            target=batch_processing,
            args=(current_app._get_current_object(), batch_text, state.batch_status),
        )
        state.batch_process.start()
    return render_template("batch.html", ready=state.batch_status.value)


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
    return jsonify({"ready": state.batch_status.value})


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
        state.batch_process.terminate()
        state.batch_process.join()
    state.batch_process = None
    state.reset_vals()
    return jsonify({"ready": state.batch_status.value})
