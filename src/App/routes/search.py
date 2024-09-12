# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

from flask import (
    Blueprint,
    render_template,
    request,
    current_app,
    jsonify,
)
import io, sys, shutil
from Comparison.Comparison import comparisonFunction as cf
import multiprocessing

search_bp = Blueprint("search", __name__)


def search_processing(
    app,
    query,
    smarts,
    smarts_num,
    finnum,
    tarray,
    incEle,
    search_status,
    captured_output,
    result_queue,
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

    This function updates the search status to "true" after processing.
    """
    with app.app_context():
        try:
            results = cf(query, finnum, tarray, incEle, smarts, smarts_num)
            result_queue.put(
                {
                    "queryval": query,
                    "results": results[0],
                    "query_models": results[1],
                    "subfailed": results[2],
                    "status": "true",
                }
            )
        except Exception as e:
            result_queue.put({"status": "false", "error": str(e)})
            print(f"Error during search processing: {e}")


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
      The order is:
        - Melting Point (MP)
        - Boiling Point (BP)
        - Log P (logP)
        - Vapor Pressure (VP)
        - Henry's Law Constant (Hlaw)
    - state.weights: An array of weights corresponding to the importance of each property in tarray.
    - state.params: Parameters for the comparison function.
    - state.results: The results of the comparison function.
    - state.query_models: Models used in the query.
    - state.subfailed: A flag indicating if substructure searching failed.
    - state.captured_output: Captured output to display progress.

    Returns:
    - GET: Renders the 'search.html' template.
    - POST: Renders the 'search.html' template with the status.
    """
    state = current_app.state
    state.reset_vals()

    if request.method == "POST":
        # Retrieve form data
        query = request.form.get("query")
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
        ]

        if query:
            # Update the status to "loading"
            state.search_status.value = "loading"
            # Terminate any existing process
            if state.search_process and state.search_process.is_alive():
                state.search_process.terminate()
                state.search_process.join()
            # Create a queue to receive results from the worker process
            result_queue = multiprocessing.Queue()
            # Start the search processing in a separate process
            state.search_process = multiprocessing.Process(
                target=search_processing,
                args=(
                    current_app._get_current_object(),
                    query,
                    smarts,
                    smarts_num,
                    state.finnum,
                    state.tarray,
                    incEle,
                    state.search_status,
                    state.captured_output,
                    result_queue,
                ),
            )
            state.search_process.start()
            # Store the result queue in the state for later access
            state.result_queue = result_queue

            # Return the search template with loading status
            return render_template("search.html", ready=state.search_status.value)

    # Render the search.html template for GET requests
    return render_template("search.html", ready=state.search_status.value)


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
    if hasattr(state, "result_queue") and not state.result_queue.empty():
        result = state.result_queue.get()
        if result["status"] == "true":
            state.queryval = result["queryval"]
            state.results = result["results"]
            state.query_models = result["query_models"]
            state.subfailed = result["subfailed"]
            state.search_status.value = "true"
        else:
            state.reset_vals()
            print(f"Error during search processing: {result.get('error')}")
        state.search_process = None

    return jsonify({"ready": state.search_status.value})


@search_bp.route("/search_cancel", methods=["POST"])
def search_cancel():
    """
    Handles the cancellation of search processing.

    This endpoint cancels the current search processing and resets the state.

    Returns:
    - JSON object with the 'ready' status set to "false".
    """
    state = current_app.state
    if state.search_process and state.search_process.is_alive():
        state.search_process.terminate()
        state.search_process.join()
    state.search_process = None
    state.reset_vals()
    return jsonify({"ready": state.search_status.value})
