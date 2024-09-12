# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pubchempy as pcp
from flask import Blueprint, jsonify, Response, current_app
from Comparison.utils.gen import parseQuery as pq

api_bp = Blueprint("api", __name__)


@api_bp.route("/api/chemicals")
def getChemicals():
    """
    Fetches chemical information for a list of CIDs stored in the application's state.

    Inputs:
    - None (uses current_app.state.cidarr for input data)

    Returns:
    - JSON response containing a list of dictionaries, each with the following keys:
        - 'cid': The chemical identifier (CID)
        - 'name': The IUPAC name or synonym of the chemical
        - 'im': URL to an image of the chemical structure
    """
    state = current_app.state
    finalres = []

    # Fetch properties for the list of CIDs
    properties = pcp.get_properties(["IUPACName"], state.cidarr)

    for i in range(len(state.cidarr)):
        property = properties[i]
        cid = state.cidarr[i]
        iupac_name = property.get("IUPACName", "")
        name = iupac_name
        if name is None:
            try:
                synonyms = pcp.get_synonyms(cid)[0]["Synonym"]
                if synonyms:
                    name = synonyms[0]
                else:
                    name = "CID: " + str(cid)
            except:
                name = "CID: " + str(cid)
        imlink = (
            "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=" + cid + "&t=l"
        )
        finalres.append({"cid": cid, "name": name, "im": imlink})
    return jsonify(finalres)


@api_bp.route("/api/query")
def query():
    """
    Parses a query value from the application's state and fetches chemical information.

    Inputs:
    - None (uses current_app.state.queryval for input data)

    Returns:
    - JSON response containing a list of dictionaries, each with the following keys:
        - 'cid': The chemical identifier (CID)
        - 'name': The name of the chemical
        - 'im': URL to an image of the chemical structure
    - If no valid result is found, returns JSON response with None.
    """
    state = current_app.state
    finalres = []
    res = pq(state.queryval)
    name = res[2]
    if res[2] != -1:
        imlink = (
            "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid="
            + str(res[1])
            + "&t=l"
        )
        finalres.append({"cid": str(res[1]), "name": name, "im": imlink})
        return jsonify(finalres)
    return jsonify(None)


@api_bp.route("/fetch_log")
def fetch_log():
    try:
        with open("logs/comparison-root.log", "r") as f:
            log_content = f.read()
        return log_content, 200
    except Exception as e:
        return str(e), 500
