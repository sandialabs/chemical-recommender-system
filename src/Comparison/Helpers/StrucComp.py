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


def extraStrucComp(comp_dict):
    """
    Performs additional structural comparisons and updates the comparison dictionary.

    This function reads thermophysical and structural data from a CSV file, performs comparisons based on
    molecular weight and other structural properties, and updates the comparison dictionary with the results.

    Inputs:
    - comp_dict: A dictionary where keys are chemical identifiers (CIDs) and values are lists containing
      comparison metrics. The dictionary is updated in place.

    Returns:
    - comp_dict: The updated comparison dictionary with adjusted structural comparison values.

    CSV File:
    - The function reads data from "src/Comparison/LocalIO/Thermout.csv" which contains the following columns:
        - "MoleculeID": The chemical identifier (CID)
        - "MolWeight": The molecular weight of the molecule
        - "nbRing": Number of rings in the molecule
        - "nbLipinskiFailures": Number of Lipinski rule violations
        - "TopoPolSurfAir": Topological polar surface area
        - "nbC": Number of carbon atoms

    Processing Steps:
    1. Read the CSV file and extract the CIDs and molecular weights.
    2. Initialize the final comparison array and set the initial comparison value for the query CID.
    3. Perform molecular weight comparisons and update the comparison dictionary.
    4. Perform additional structural comparisons based on predefined weights and properties.
    5. Update the comparison dictionary with the adjusted structural comparison values.
    """
    logger.info("Starting extra structural comparisons.")

    try:
        # Read the CSV file
        df = pd.read_csv("src/Comparison/LocalIO/Thermout.csv")
        cids = df["MoleculeID"].values
        finalarr = [1] * len(cids)
        qcid = cids[0]
        comp_dict[str(qcid)][2] = 1
        cids = cids[1:]

        # Perform molecular weight comparison
        vals = df["MolWeight"].values
        comp = vals[0]
        if float(comp) == 0:
            comp = float(comp) + 0.001
        vals = vals[1:]
        for i in range(len(cids)):
            key = str(cids[i])
            finalval = 1 / (1 + abs((float(vals[i]) - comp) / (comp)))
            comp_dict[key][2] = finalval

        # Define weights and structural properties for additional comparisons
        weights = [0.1, 0.1, 0.1]
        comparisons = ["nbRing", "nbLipinskiFailures", "TopoPolSurfAir", "nbC"]

        for j in range(len(weights)):
            # Compare structural properties
            vals = df[comparisons[j]].values
            comp = vals[0]
            vals = vals[1:]
            # Adjust values of the finalarr array with coefficients for scaling on thermal comparison
            for i in range(len(cids)):
                if float(comp) == 0:
                    finalarr[i] *= (1 / (1 + abs((float(vals[i]) / 0.001)))) ** weights[
                        j
                    ]
                else:
                    finalarr[i] *= (
                        1 / (1 + abs((float(vals[i]) - comp) / comp))
                    ) ** weights[j]

            # Adjust structural comparison values and update the aggregate value
            for i in range(len(cids)):
                key = str(cids[i])
                comp_dict[key][1] *= finalarr[i]

        logger.info("Extra structural comparisons completed successfully.")
        return comp_dict
    except Exception as e:
        logger.error(f"Error during extra structural comparisons: {e}")
        raise
