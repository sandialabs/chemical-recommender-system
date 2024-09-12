# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pandas as pd
import pubchempy as pcp


def compName(cid):
    """
    Get the compound name for a given CID, with a fallback to the CID if no name is found.

    Inputs:
    - cid: The chemical identifier (CID).

    Returns:
    - val: A string containing the CID and the compound name (truncated if necessary).
    """
    try:
        val = "CID: " + str(cid)
        properties = pcp.get_properties(["IUPACName"], cid)[0]
        name = properties.get("IUPACName", None)

        if name is None:
            synonyms = pcp.get_synonyms(cid)[0].get("Synonym", [])
            name = synonyms[0] if synonyms else "CID: " + str(cid)

        if len(name) > 20:
            val = val + ": " + name[:18] + "..."
        else:
            val = val + ": " + name

        return val
    except:
        return "CID: " + str(cid)


def createCSVOut(results, tarray, containers=[], query_models=[]):
    """
    Create a CSV output file from the results and user search parameters.

    Inputs:
    - results: A list of dictionaries containing the results data.
    - tarray: An array of boolean values indicating which thermophysical properties are to be included.
      The order is:
        - Melting Point (MP)
        - Boiling Point (BP)
        - Log P (logP)
        - Vapor Pressure (VP)
        - Henry's Law Constant (Hlaw)
    - containers: A list of container names for additional columns (default is an empty list).
    - query_models: A list of query models corresponding to the containers (default is an empty list).

    Returns:
    - None (saves the CSV file to "src/App/static/LocalIO/results.csv")
    """
    df = pd.DataFrame(data=results, index=None)

    # Initialize columns for CSV based on user's search parameters
    columns = [
        "Molecule Name",
        "Final Similarity Score",
        "Molecular Weight Similarity",
        "Structural Similarity",
        "Thermophysical Similarity",
        "Predicted Toxicity",
        "SA Score",
    ]
    if containers:
        for i in range(len(containers)):
            if query_models[i] is not None:
                columns.append(containers[i])
    if tarray[0]:
        columns.append("MP_pred")
    if tarray[1]:
        columns.append("BP_pred")
    if tarray[2]:
        columns.append("LogP_pred")
    if tarray[3]:
        columns.append("LogHL_pred")
    if tarray[4]:
        columns.append("LogVP_pred")
    columns.append("Predicted BCF")
    columns.append("CATMoS EPA")
    columns.append("LD50")
    df.columns = columns

    # Transform molecule name column from CIDs to actual names
    try:
        df["Molecule Name"] = df["Molecule Name"].apply(compName)
    except:
        None

    # Specify the query and round all values to 3 decimals
    pd.options.mode.chained_assignment = None  # default='warn'
    df["Molecule Name"].iloc[0] = "QUERY: " + df["Molecule Name"].iloc[0]
    pd.options.mode.chained_assignment = "warn"  # Reset to default

    df = df.round(3)

    # Save CSV into LocalIO
    df.to_csv("src/App/static/LocalIO/results.csv", index=False)
