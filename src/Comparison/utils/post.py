# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pandas as pd

from Comparison.utils.norm import normalizeDict as nd
from Comparison.utils.csv import createCSVOut as co


def createCheckstring(tarray):
    """
    Create a list of prediction columns based on the input array.

    Inputs:
    - tarray: An array of boolean values indicating which thermophysical properties are to be included.
      The order is:
        - Melting Point (MP)
        - Boiling Point (BP)
        - Log P (logP)
        - Vapor Pressure (VP)
        - Henry's Law Constant (Hlaw)

    Returns:
    - checkstring: A list of prediction column names.
    """
    checkstring = []
    if tarray[0]:
        checkstring.append("MP_pred")  # Melting Point prediction
    if tarray[1]:
        checkstring.append("BP_pred")  # Boiling Point prediction
    if tarray[2]:
        checkstring.append("LogP_pred")  # LogP prediction
    if tarray[3]:
        checkstring.append("LogHL_pred")  # Henry's Law prediction
    if tarray[4]:
        checkstring.append("LogVP_pred")  # Vapor Pressure prediction
    checkstring.append("LogBCF_pred")  # Bioconcentration Factor prediction
    checkstring.append("CATMoS_EPA_pred")  # EPA prediction
    checkstring.append("CATMoS_LD50_pred")  # LD50 prediction
    return checkstring


def processDict(results, weights, qcid, finnum, query_models=None, containers=None):
    """
    Process and normalize the results dictionary.

    Inputs:
    - results: A dictionary containing the results data.
    - weights: A list of weights for the query.
    - qcid: The query chemical identifier (CID).
    - finnum: The number of final results to be returned.
    - query_models: A list of query models (default is None).
    - containers: A list of container names for additional columns (default is None).

    Returns:
    - final_results: The processed and normalized results dictionary.
    """
    num_containers = 0
    if containers:
        num_containers = len(containers)

    if query_models:
        for val in query_models:
            if val is None:
                num_containers -= 1

    # Set the comparison values to None for the query on itself
    results[str(qcid)][0] = None
    results[str(qcid)][1] = None
    results[str(qcid)][2] = None
    results[str(qcid)][3] = None

    # Normalize the results dictionary
    results = nd(results)

    # Calculate the combined score for each compound
    for key, val in results.items():
        if val[0] is None:
            None
        else:
            val[0] = (
                (val[1] ** float(weights[0]))  # Structural similarity weight
                * (
                    val[2] ** (float(weights[1]) / 10)
                )  # Molecular weight similarity weight
                * (val[3] ** float(weights[2]))  # Thermophysical similarity weight
                * (val[5] ** float(weights[4]))  # Synthetic availability weight
                / (val[4] ** float(weights[3]))  # Toxicity weight
            )
            for i in range(num_containers):
                val[0] *= val[i + 5] ** float(
                    weights[i + 5]
                )  # Additional weights for containers

    # Hold the query compound's key and value
    hold_key = str(qcid)
    hold_val = results[str(qcid)]
    results.pop(hold_key)

    # Sort the results dictionary by the combined score in descending order
    results = dict(sorted(results.items(), key=lambda x: x[1][0], reverse=True))
    results = {
        k: results[k] for k in list(results)[:finnum]
    }  # Keep only the top 'finnum' results

    # Add the query compound back to the final results
    final_results = {hold_key: hold_val}
    for key, val in results.items():
        final_results[key] = val

    return final_results


def produceCSV(processed_dict, checkstring, tarray, containers=None, query_models=None):
    """
    Produce a CSV file from the processed dictionary.

    Inputs:
    - processed_dict: The processed results dictionary.
    - checkstring: A list of prediction column names.
    - tarray: An array of boolean values indicating which thermophysical properties are to be included.
    - containers: A list of container names for additional columns (default is None).
    - query_models: A list of query models (default is None).

    Returns:
    - finalresults: A list of lists containing the final results data.
    - fullcids: A list of all compound IDs.
    """
    num_containers = 0
    if containers:
        num_containers = len(containers)

    # Read the CSV file into a DataFrame
    df = pd.read_csv("src/Comparison/LocalIO/Thermout.csv")
    finalresults = []

    # Iterate over the processed dictionary to create the final results list
    for key, val in processed_dict.items():
        newarr = []
        newarr.append(key)  # Add the compound ID
        for i in range(6 + num_containers):
            newarr.append(val[i])  # Add the normalized values

        cid = int(key)
        matching_row = df.index[df["MoleculeID"] == cid].tolist()[0]
        for colname in checkstring:
            newarr.append(
                str(df.loc[matching_row, colname])
            )  # Add the prediction values

        finalresults.append(newarr)

    fullcids = []
    for i in range(len(finalresults)):
        fullcids.append(finalresults[i][0])  # Collect all compound IDs

    # Save finalresults into the correct CSV
    co(finalresults, tarray, containers, query_models)
    return finalresults, fullcids


def formatResults(finalresults, qcid):
    """
    Format the final results for display.

    Inputs:
    - finalresults: A list of lists containing the final results data.
    - qcid: The query chemical identifier (CID).

    Returns:
    - finalresults: The formatted final results list.
    - cidarr: A list of unique compound IDs.
    """
    finalresults = finalresults[:21]  # Keep only the top 21 results to show in GUI
    cidarr = []
    for i in range(len(finalresults)):
        if finalresults[i][0] not in cidarr and finalresults[i][0] != str(qcid):
            cidarr.append(finalresults[i][0])  # Collect unique compound IDs
        if i == 0:
            finalresults[i][0] = (
                "<strong>" + "Query:" + "</strong>"
            )  # Mark the query compound
        else:
            finalresults[i][0] = (
                "<strong>"
                + str(i)
                + "</strong>"
                + ". "
                + finalresults[i][0]  # Format the result index and ID
            )
    return finalresults, cidarr
