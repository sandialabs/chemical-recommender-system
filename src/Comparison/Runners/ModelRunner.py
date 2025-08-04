# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import requests

# Import centralized logging
from utils.progress_logger import get_progress_logger

# Get a logger for this module (will be replaced with job-specific logger)
logger = None

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


def call_service(smiles_list, image, logger=None):
    """
    Call an external service to compute results based on a list of SMILES strings.

    Inputs:
    - smiles_list: A list of SMILES strings to be processed.
    - image: The hostname or IP address of the service.
    - logger: Logger instance for debugging (optional).

    Returns:
    - The JSON response from the service.
    """
    url = "http://" + str(image) + ":3012/compute"
    payload = {"smiles": smiles_list}

    if logger:
        logger.info(f"Attempting to call model service at: {url}")
        logger.info(f"Service payload contains {len(smiles_list)} SMILES strings")
        logger.debug(f"Full payload: {payload}")
    else:
        print(f"WARNING: No logger provided to call_service - calling {url}")

    try:
        # Explicitly set proxies to None to ignore any proxy settings
        if logger:
            logger.info("Sending POST request to model service...")
        
        response = requests.post(url, json=payload, proxies={"http": None, "https": None})
        
        if logger:
            logger.info(f"Model service responded with status code: {response.status_code}")
            
        response.raise_for_status()  # Raise an exception for bad status codes
        
        response_json = response.json()
        if logger:
            logger.info("Successfully received and parsed JSON response from model service")
            logger.debug(f"Response keys: {list(response_json.keys())}")
            logger.debug(f"Full service response: {response_json}")
            
        return response_json
        
    except requests.exceptions.ConnectionError as e:
        error_msg = f"Connection error when calling model service {url}: {str(e)}"
        if logger:
            logger.error(error_msg)
        else:
            print(f"ERROR: {error_msg}")
        raise
    except requests.exceptions.HTTPError as e:
        error_msg = f"HTTP error from model service {url}: {response.status_code} - {str(e)}"
        if logger:
            logger.error(error_msg)
        else:
            print(f"ERROR: {error_msg}")
        raise
    except requests.exceptions.RequestException as e:
        error_msg = f"Request error when calling model service {url}: {str(e)}"
        if logger:
            logger.error(error_msg)
        else:
            print(f"ERROR: {error_msg}")
        raise
    except ValueError as e:
        error_msg = f"JSON parsing error from model service response: {str(e)}"
        if logger:
            logger.error(error_msg)
            logger.error(f"Raw response content: {response.text}")
        else:
            print(f"ERROR: {error_msg}")
        raise
    except Exception as e:
        error_msg = f"Unexpected error when calling model service {url}: {str(e)}"
        if logger:
            logger.error(error_msg)
        else:
            print(f"ERROR: {error_msg}")
        raise


def runModels(comp_dict, image, query_smi, smiles_dict, job_id="default"):
    """
    Run models on a list of SMILES strings and update the comparison dictionary.

    Inputs:
    - comp_dict: A dictionary where keys are chemical identifiers (CIDs) and values are lists containing comparison metrics.
    - image: The hostname or IP address of the service.
    - query_smi: The query SMILES string.
    - smiles_dict: A dictionary mapping CIDs to SMILES strings.
    - job_id: Unique identifier for this job (for logging purposes).

    Returns:
    - A list containing the updated comparison dictionary and the computed result for the query SMILES string.
    """
    logger = get_progress_logger(job_id)
    logger.info("============ Starting Model Runner ============")
    logger.info(f"Running models on SMILES strings using service: {image}")
    logger.info(f"Query SMILES: {query_smi}")
    logger.info(f"Number of compounds to process: {len(comp_dict)}")
    logger.info(f"Available SMILES in smiles_dict: {len(smiles_dict)}")
    
    key_arr = [query_smi]
    keys = []
    
    logger.info("Building SMILES array for model service...")
    for key, val in comp_dict.items():
        keys.append(key)
        if key == "-1":
            logger.debug(f"Using query SMILES for key {key}")
            key_arr.append(query_smi)
        else:
            if key in smiles_dict:
                key_arr.append(smiles_dict[key])
                logger.debug(f"Added SMILES for CID {key}: {smiles_dict[key]}")
            else:
                logger.warning(f"CID {key} not found in smiles_dict")
                key_arr.append("")  # Add empty string as placeholder

    logger.info(f"Final SMILES array contains {len(key_arr)} entries")
    
    try:
        result = call_service(key_arr, image, logger)
        logger.info("Model service call completed successfully")
    except Exception as e:
        logger.error(f"Model service call failed: {str(e)}")
        raise

    # Validate the result
    if not isinstance(result, dict):
        logger.error(f"Model service returned unexpected type: {type(result)}, expected dict")
        return [comp_dict, None]
        
    logger.info(f"Model service returned results for {len(result)} SMILES")
    logger.debug(f"Result keys: {list(result.keys())}")

    if None in result.values():
        logger.warning("Model service returned None values for some candidates")
        none_keys = [k for k, v in result.items() if v is None]
        logger.warning(f"SMILES with None results: {none_keys}")
        logger.warning("Added Model does not return a value for all candidates, ignoring model")
        return [comp_dict, None]

    # Check if query SMILES is in results
    if query_smi not in result:
        logger.error(f"Query SMILES '{query_smi}' not found in model results")
        logger.error(f"Available keys in result: {list(result.keys())}")
        return [comp_dict, None]

    # Use the output array in the rest of your Python code
    finalarr = [1] * len(keys)
    comp = result[query_smi]
    logger.info(f"Query compound model result: {comp}")

    # Adjust values of the finalarr array with coefficients for scaling on thermal comparison
    logger.info("Processing model results for candidate compounds...")
    vals = []
    missing_smiles = []
    
    for i, key in enumerate(keys):
        if key == "-1":
            smiles = query_smi
        else:
            smiles = smiles_dict.get(key, "")
            
        if smiles in result:
            converted_val = convertItem(result[smiles])
            vals.append(converted_val)
            logger.debug(f"CID {key}: SMILES '{smiles}' -> result {result[smiles]} -> converted {converted_val}")
        else:
            logger.warning(f"CID {key} with SMILES '{smiles}' not found in model results")
            missing_smiles.append((key, smiles))
            vals.append(None)
    
    if missing_smiles:
        logger.warning(f"Missing results for {len(missing_smiles)} compounds: {missing_smiles}")
    
    comp = convertItem(comp)
    logger.info(f"Converted query compound result: {comp}")

    if comp == 0:
        logger.warning("Query compound result is 0, adjusting to 0.001 to avoid division by zero")
        comp += 0.001

    logger.info("Calculating similarity scores based on model results...")
    valid_results = 0
    
    for i in range(len(keys)):
        if vals[i] is not None:
            try:
                similarity_score = 1 / (1 + abs((vals[i] - comp) / comp))
                finalarr[i] *= similarity_score
                valid_results += 1
                logger.debug(f"CID {keys[i]}: val={vals[i]}, comp={comp}, similarity={similarity_score}")
            except (TypeError, ZeroDivisionError) as e:
                logger.error(f"Error calculating similarity for CID {keys[i]}: {e}")
                finalarr[i] = 0  # Set to 0 if calculation fails
        else:
            logger.warning(f"CID {keys[i]} has None value, setting similarity to 0")
            finalarr[i] = 0

    logger.info(f"Successfully calculated similarities for {valid_results}/{len(keys)} compounds")

    logger.info("Updating comparison dictionary with model results...")
    updated_count = 0
    
    for i in range(len(keys)):
        key_str = str(keys[i])
        if key_str not in comp_dict:
            logger.error(f"Key {key_str} not found in comp_dict during model update")
            continue
        if len(comp_dict[key_str]) < 1:
            logger.error(f"comp_dict[{key_str}] has insufficient elements: {comp_dict[key_str]}")
            continue
            
        original_score = comp_dict[key_str][0]
        comp_dict[key_str][0] *= finalarr[i]
        comp_dict[key_str].append(finalarr[i])
        updated_count += 1
        
        logger.debug(f"Updated CID {key_str}: original_score={original_score}, "
                    f"model_factor={finalarr[i]}, new_score={comp_dict[key_str][0]}")

    logger.info(f"Successfully updated {updated_count}/{len(keys)} compounds in comparison dictionary")
    logger.info("============ Model Runner Completed ============")
    logger.debug(f"Final comparison dictionary sample: {dict(list(comp_dict.items())[:3])}")
    
    return [comp_dict, comp]
