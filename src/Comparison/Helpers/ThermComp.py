# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pandas as pd
import logging

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


def thermalComparison(comp_dict, tarray):
    """
    Performs thermal property comparisons and updates the comparison dictionary.

    This function reads thermophysical data from a CSV file, performs comparisons based on user-selected
    thermophysical properties, and updates the comparison dictionary with the results.

    Inputs:
    - comp_dict: A dictionary where keys are chemical identifiers (CIDs) and values are lists containing
      comparison metrics. The dictionary is updated in place.
    - tarray: An array of boolean values indicating which thermophysical properties are to be calculated.
      The order is:
        - Melting Point (MP)
        - Boiling Point (BP)
        - Log P (logP)
        - Vapor Pressure (VP)
        - Henry's Law Constant (Hlaw)

    Returns:
    - comp_dict: The updated comparison dictionary with adjusted thermal comparison values.

    CSV File:
    - The function reads data from "src/Comparison/LocalIO/Thermout.csv" which contains the following columns:
        - "MoleculeID": The chemical identifier (CID)
        - "MP_pred": Predicted melting point
        - "BP_pred": Predicted boiling point
        - "LogP_pred": Predicted log P
        - "LogHL_pred": Predicted Henry's Law constant
        - "LogVP_pred": Predicted vapor pressure

    Processing Steps:
    1. Read the CSV file and extract the CIDs.
    2. Initialize the final comparison array.
    3. Define the endpoints for thermophysical properties based on user input.
    4. Perform comparisons for each selected endpoint and update the final comparison array.
    5. Update the comparison dictionary with the adjusted thermal comparison values.
    """
    logger.info("Starting thermal property comparisons.")

    try:
        # Define which thermophysical comparisons to do based on user input
        MP = tarray[0]
        BP = tarray[1]
        logP = tarray[2]
        Hlaw = tarray[3]
        VP = tarray[4]

        # The output goes to Thermout and this is read in using pandas and stored in df
        df = pd.read_csv("src/Comparison/LocalIO/Thermout.csv")
        cids = df["MoleculeID"].values
        finalarr = [1] * len(cids)
        cids = cids[1:]

        # Define the endpoints for thermophysical properties
        endpoints = [
            [MP, "MP_pred"],
            [BP, "BP_pred"],
            [logP, "LogP_pred"],
            [Hlaw, "LogHL_pred"],
            [VP, "LogVP_pred"],
        ]

        # If each of the endpoints is on, they are considered in computation
        for item in endpoints:
            if item[0]:
                vals = df[item[1]].values
                comp = vals[0]
                if float(comp) == 0:
                    comp = float(comp) + 0.001
                vals = vals[1:]
                # Adjust values of the finalarr array with coefficients for scaling on thermal comparison
                for i in range(len(cids)):
                    finalarr[i] *= 1 / (1 + abs((float(vals[i]) - comp) / (comp)))

        # Add to the end of the array and adjust the aggregate value
        for i in range(len(cids)):
            comp_dict[str(cids[i])][3] = finalarr[i]

        logger.info("Thermal property comparisons completed successfully.")
        return comp_dict
    except Exception as e:
        logger.error(f"Error during thermal property comparisons: {e}")
        raise
