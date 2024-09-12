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


def toxicComparison(comp_dict):
    """
    Performs toxicity property comparisons and updates the comparison dictionary.

    This function reads toxicity data from a CSV file, performs calculations based on specific toxicity metrics,
    and updates the comparison dictionary with the results.

    Inputs:
    - comp_dict: A dictionary where keys are chemical identifiers (CIDs) and values are lists containing
      comparison metrics. The dictionary is updated in place.

    Returns:
    - comp_dict: The updated comparison dictionary with adjusted toxicity comparison values.

    CSV File:
    - The function reads data from "src/Comparison/LocalIO/Thermout.csv" which contains the following columns:
        - "MoleculeID": The chemical identifier (CID)
        - "LogBCF_pred": Predicted bioconcentration factor (BCF)
        - "CATMoS_EPA_pred": Predicted EPA category
        - "CATMoS_LD50_pred": Predicted LD50 value

    Processing Steps:
    1. Read the CSV file and extract the relevant toxicity data.
    2. For each CID in the comparison dictionary, find the matching row in the CSV file.
    3. Extract the BCF, EPA, and LD50 values for the CID.
    4. Calculate the total toxicity value using the formula: totalval = BCF * EPA * LD50 / 2.
    5. Update the comparison dictionary with the calculated toxicity value.
    """
    logger.info("Starting toxicity property comparisons.")

    try:
        # Read the CSV file
        df = pd.read_csv("src/Comparison/LocalIO/Thermout.csv")

        # Iterate over each CID in the comparison dictionary
        for cid, vals in comp_dict.items():
            cid = int(cid)
            matching_row = df.index[df["MoleculeID"] == cid].tolist()[0]

            # Extract toxicity values from the CSV file
            BCF = float(df.loc[matching_row, "LogBCF_pred"])
            EPA = 1 + float(df.loc[matching_row, "CATMoS_EPA_pred"])
            LD50 = float(df.loc[matching_row, "CATMoS_LD50_pred"]) / 1000

            # Calculate the total toxicity value
            totalval = BCF * EPA * LD50 / 2

            # Update the comparison dictionary with the calculated toxicity value
            vals[4] = totalval

        logger.info("Toxicity property comparisons completed successfully.")
        return comp_dict
    except Exception as e:
        logger.error(f"Error during toxicity property comparisons: {e}")
        raise
