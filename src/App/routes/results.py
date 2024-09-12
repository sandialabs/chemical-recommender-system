# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

from flask import Blueprint, render_template, request, current_app
import pandas, io, sys, shutil
from Comparison.utils.gen import parseQuery as pq
from Comparison.utils.graph import graphResults as gr
from Comparison.utils.pdf import combinePDFs as cp
from Comparison.utils.post import (
    createCheckstring as cc,
    processDict as pd,
    produceCSV as pc,
    formatResults as fr,
)

results_bp = Blueprint("results", __name__)


@results_bp.route("/results", methods=["GET", "POST"])
def resultsFunc():
    """
    Processes and displays the results based on user input and application state.

    This function handles both GET and POST requests to the /results route. It processes user-provided weights,
    generates results, and renders the results.html template with the final results.

    Inputs:
    - None directly (uses current_app.state and request.form for input data)

    Application State:
    - state.tarray: An array of boolean values indicating which thermophysical properties are to be calculated.
      The order is:
        - Melting Point
        - Boiling Point
        - Log P
        - Vapor Pressure
        - Henry's Law Constant
    - state.weights: An array of weights corresponding to the importance of each property in tarray.
    - state.queryval: The query value provided by the user.
    - state.results: The initial results to be processed.
    - state.finnum: A parameter used in processing the results.
    - state.params: Additional parameters for generating the final report.
    - state.subfailed: A flag indicating if any subprocesses failed.

    Returns:
    - Renders the 'results.html' template with the following context variables:
        - 'finalresults': The final processed results to be displayed
        - 'tarray': The array of target values from the application state
        - 'weights': The array of weights from the application state
    """
    state = current_app.state

    # Handle POST request to update weights
    if request and request.method == "POST":
        state.weights[0] = request.form.get("weight1")
        state.weights[1] = request.form.get("weight2")
        state.weights[2] = request.form.get("weight3")
        state.weights[3] = request.form.get("weight4")
        state.weights[4] = request.form.get("weight5")

    # Create a checkstring based on the target array
    checkstring = cc(state.tarray)
    # Parse the query value to get the query CID
    qcid = (pq(state.queryval))[1]
    if qcid is None:
        qcid = -1

    # Process the results dictionary with the given weights and query CID
    processed_dict = pd(state.results, state.weights, qcid, state.finnum)

    # Save the processed dictionary to a CSV file
    dfcsv = pandas.DataFrame(processed_dict)
    dfcsv.to_csv(r"src/Comparison/LocalIO/data.csv", index=False)

    shutil.copy(
        "src/Comparison/LocalIO/Thermout.csv",
        "src/App/static/LocalIO/OPERA-Predictions.csv",
    )

    # Generate graphs based on the query CID or query value
    if qcid != -1:
        gr(qcid)
    else:
        gr(qcid, state.queryval)

    # Produce the final results and full CIDs list
    finalresults, fullcids = pc(
        processed_dict, checkstring, state.tarray, containers=None, query_models=None
    )

    # Combine PDFs for the final report
    cp(
        fullcids,
        state.queryval,
        qcid,
        state.params,
        state.weights,
        True,
        state.subfailed,
    )

    # Format the final results and update the CID array in the state
    finalresults, state.cidarr = fr(finalresults, qcid)

    # Render the results.html template with the final results and state variables
    return render_template(
        "results.html",
        finalresults=finalresults,
        tarray=state.tarray,
        weights=state.weights,
    )
