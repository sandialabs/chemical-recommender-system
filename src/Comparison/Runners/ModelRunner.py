# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import requests
import os
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


def convertItem(item):
    """
    Convert an item to a float if possible, otherwise handle it appropriately.

    Inputs:
    - item: The item to be converted.

    Returns:
    - The converted item as a float, None if the item is "None", or the item itself if conversion fails.
    """
    try:
        # Try converting to float
        return float(item)
    except ValueError:
        # If conversion fails, check if the item is 'None'
        if item == "None":
            return None
        else:
            # If it's not 'None', return the item as is or handle differently
            return item


def call_service(smiles_list, image):
    """
    Call an external service to compute results based on a list of SMILES strings.

    Inputs:
    - smiles_list: A list of SMILES strings to be processed.
    - image: The hostname or IP address of the service.

    Returns:
    - The JSON response from the service.
    """
    url = "http://" + str(image) + ":3012/compute"
    payload = {"smiles": smiles_list}

    logger.debug(f"Calling service at {url} with payload: {payload}")

    # Explicitly set proxies to None to ignore any proxy settings
    response = requests.post(url, json=payload, proxies={"http": None, "https": None})

    logger.debug(f"Service response: {response.json()}")
    return response.json()


def runModels(comp_dict, image, query_smi, smiles_dict):
    """
    Run models on a list of SMILES strings and update the comparison dictionary.

    Inputs:
    - comp_dict: A dictionary where keys are chemical identifiers (CIDs) and values are lists containing comparison metrics.
    - image: The hostname or IP address of the service.
    - query_smi: The query SMILES string.
    - smiles_dict: A dictionary mapping CIDs to SMILES strings.

    Returns:
    - A list containing the updated comparison dictionary and the computed result for the query SMILES string.
    """
    logger.info("Running models on SMILES strings.")
    key_arr = [query_smi]
    keys = []
    for key, val in comp_dict.items():
        keys.append(key)
        if key == "-1":
            key_arr.append(query_smi)
        else:
            key_arr.append(smiles_dict[key])

    result = call_service(key_arr, image)

    if None in result.values():
        logger.warning(
            "Added Model does not return a value for all candidates, ignoring model"
        )
        return [comp_dict, None]

    # Use the output array in the rest of your Python code
    finalarr = [1] * len(keys)
    comp = result[query_smi]

    # Adjust values of the finalarr array with coefficients for scaling on thermal comparison
    vals = [convertItem(result[smiles_dict[key]]) for key in keys]
    comp = convertItem(comp)

    if comp == 0:
        comp += 0.001

    for i in range(len(keys)):
        finalarr[i] *= 1 / (1 + abs((vals[i] - comp) / comp))

    for i in range(len(keys)):
        comp_dict[str(keys[i])][0] *= finalarr[i]
        comp_dict[str(keys[i])].append(finalarr[i])

    logger.debug(f"Updated comparison dictionary: {comp_dict}")
    return [comp_dict, comp]
