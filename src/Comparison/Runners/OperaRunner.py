# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pubchempy as pcp
import subprocess
import threading
import logging

# Ensure logging is configured in the main script
logger = logging.getLogger(__name__)

# Add handlers to the logger to ensure logs go to both comparison.log and comparison-root.log
comparison_handler = logging.FileHandler("logs/comparison.log", mode="a")
comparison_handler.setLevel(logging.DEBUG)
comparison_formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
comparison_handler.setFormatter(comparison_formatter)

root_handler = logging.FileHandler("logs/comparison-root.log", mode="a")
root_handler.setLevel(logging.INFO)
root_formatter = logging.Formatter("%(message)s")
root_handler.setFormatter(root_formatter)

logger.addHandler(comparison_handler)
logger.addHandler(root_handler)


def read_output(process):
    """
    Capture OPERA output to stdout to display in GUI.

    Inputs:
    - process: The subprocess.Popen object representing the running process.

    Returns:
    - None
    """
    for line in iter(process.stdout.readline, b""):
        stripped_line = line.strip()
        if stripped_line:  # Only log non-empty lines
            logger.info(stripped_line)
        if "molecules predicted. Total process time" in stripped_line:
            return


def run_command(command):
    """
    Use the subprocess module to send a command line like run for OPERA.

    Inputs:
    - command: The command to be executed as a string.

    Returns:
    - None
    """
    logger.info(f"Running OPERA")
    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, text=True, universal_newlines=True
    )
    output_thread = threading.Thread(target=read_output, args=(process,))
    output_thread.start()

    process.wait()
    output_thread.join()


def runOpera(comp_dict, querycid, tarray, query_smi, smiles_dict):
    """
    Run the OPERA tool to predict thermophysical properties based on user input.

    Inputs:
    - comp_dict: A dictionary where keys are chemical identifiers (CIDs) and values are lists containing comparison metrics.
    - querycid: The query chemical identifier (CID).
    - tarray: An array of boolean values indicating which thermophysical properties are to be calculated.
      The order is:
        - Melting Point (MP)
        - Boiling Point (BP)
        - Log P (logP)
        - Vapor Pressure (VP)
        - Henry's Law Constant (Hlaw)
    - query_smi: The query SMILES string.
    - smiles_dict: A dictionary mapping CIDs to SMILES strings.

    Returns:
    - None
    """
    logger.info("============ Entering OPERA Models ============")
    # Define which thermophysical comparisons to do based on user input, define endpoint string
    MP = tarray[0]
    BP = tarray[1]
    logP = tarray[2]
    Hlaw = tarray[3]
    VP = tarray[4]

    endpoints = "CATMoS BCF "
    if MP:
        endpoints += "MP"
    if BP:
        endpoints += " BP"
    if logP:
        endpoints += " logP"
    if Hlaw:
        endpoints += " HL"
    if VP:
        endpoints += " VP"
    endpoints += " StrP"

    # OPERA requires SMILES input, create a CID and SMILES array from the existing dict
    logger.info("Initializing Opera")
    if query_smi:
        querysmiles = query_smi
    else:
        querysmiles = pcp.get_properties("CanonicalSMILES", querycid)[0][
            "CanonicalSMILES"
        ]
    smilesarr = [querysmiles]
    cidarr = [querycid]
    for key, _ in comp_dict.items():
        cidarr.append(key)
        if key == "-1":
            smilesarr.append(query_smi)
        else:
            smilesarr.append(smiles_dict[key])

    # Write the SMILES and CID data to a .smi file in LocalIO
    with open("src/Comparison/LocalIO/Thermin.smi", "w") as f:
        for s, c in zip(smilesarr, cidarr):
            f.write(f"{s}\t{c}\n")

    logger.info("Starting up Opera")

    command = f"/usr/OPERA/application/OPERA --SMI src/Comparison/LocalIO/Thermin.smi -o src/Comparison/LocalIO/Thermout.csv -c -e {endpoints} -v 1"
    run_command(command)
