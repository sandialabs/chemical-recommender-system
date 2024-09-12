# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

# Add in the app.py file and connect your model as the process_smiles function

from flask import Flask, request, jsonify
from function import process_smiles

app = Flask(__name__)


@app.route("/compute", methods=["POST"])
def compute():
    data = request.json
    smiles_list = data.get("smiles")
    result = {}
    for smile in smiles_list:
        print(smile)
        try:
            result[smile] = process_smiles(smile)
        except:
            result[smile] = None
    return jsonify(result)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=3012)
