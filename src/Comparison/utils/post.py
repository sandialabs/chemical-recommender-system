# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pandas as pd
import logging
import traceback

from Comparison.utils.norm import normalizeDict as nd
from Comparison.utils.csv import createCSVOut as co

# Get logger for this module
logger = logging.getLogger(__name__)


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
    # === TECHNICAL VALIDATION ===
    if not isinstance(results, dict):
        raise ValueError(f"results must be a dictionary, got: {type(results)}")
    if not isinstance(weights, (list, tuple)):
        raise ValueError(f"weights must be a list or tuple, got: {type(weights)}")
    if not isinstance(finnum, int) or finnum <= 0:
        raise ValueError(f"finnum must be a positive integer, got: {finnum}")
    
    num_containers = 0
    if containers:
        num_containers = len(containers)

    if query_models:
        for val in query_models:
            if val is None:
                num_containers -= 1

    # Validate weights array length BEFORE using indices
    expected_weights_length = 5 + num_containers
    if len(weights) != expected_weights_length:
        raise ValueError(f"weights array length ({len(weights)}) doesn't match expected length ({expected_weights_length}). Need 5 base weights + {num_containers} container weights")
    
    # Validate that qcid exists in results
    qcid_str = str(qcid)
    if qcid_str not in results:
        raise ValueError(f"Query CID {qcid} not found in results dictionary. Available CIDs: {list(results.keys())[:5]}...")
    
    # Validate all result values have correct structure
    expected_value_length = 6 + num_containers
    for cid, values in results.items():
        if not isinstance(values, (list, tuple)):
            raise ValueError(f"Result values must be lists/tuples, CID {cid} has: {type(values)}")
        if len(values) != expected_value_length:
            raise ValueError(f"Result values must have {expected_value_length} elements, CID {cid} has {len(values)}: {values}")
    # === END VALIDATION ===

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
            # === CRITICAL FIX: Validate indices before accessing ===
            # Check that we have enough values before accessing indices
            if len(val) < 6:
                raise ValueError(f"Insufficient values in result for CID {key}: has {len(val)}, need at least 6")
            if len(val) < 6 + num_containers:
                raise ValueError(f"Insufficient values for containers in result for CID {key}: has {len(val)}, need {6 + num_containers}")
                
            # Validate weight indices before using them
            for i in range(5):
                if i >= len(weights):
                    raise IndexError(f"Weight index {i} out of bounds, weights length: {len(weights)}")
            for i in range(num_containers):
                weight_idx = i + 5
                if weight_idx >= len(weights):
                    raise IndexError(f"Container weight index {weight_idx} out of bounds, weights length: {len(weights)}")
                val_idx = i + 6
                if val_idx >= len(val):
                    raise IndexError(f"Container value index {val_idx} out of bounds, value length: {len(val)}")
            
            val[0] = (
                (val[1] ** float(weights[0]))  # Structural similarity weight
                * (
                    val[2] ** (float(weights[1]))
                )  # Molecular weight similarity weight
                * (val[3] ** float(weights[2]))  # Thermophysical similarity weight
                * (val[5] ** float(weights[4]))  # Synthetic availability weight
                / (val[4] ** float(weights[3]))  # Toxicity weight
            )
            for i in range(num_containers):
                val[0] *= val[i + 6] ** float(
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

    # Sort processed_dict by final similarity score (index 0) in descending order, excluding query (None values)
    sorted_items = []
    query_item = None
    
    for key, val in processed_dict.items():
        if val[0] is None:  # This is the query compound
            query_item = (key, val)
        else:
            sorted_items.append((key, val))
    
    # Sort by final similarity score (val[0]) in descending order
    sorted_items.sort(key=lambda x: x[1][0], reverse=True)
    
    # Add query at the beginning if it exists
    if query_item:
        sorted_items.insert(0, query_item)

    # Iterate over the sorted results to create the final results list
    for key, val in sorted_items:
        # === VALIDATION: Check val structure before accessing indices ===
        if not isinstance(val, (list, tuple)):
            logger.error(f"Invalid val type for CID {key}: {type(val)}")
            continue
        expected_val_length = 6 + num_containers
        if len(val) < expected_val_length:
            logger.error(f"Insufficient val elements for CID {key}: has {len(val)}, need {expected_val_length}")
            continue
            
        newarr = []
        # Convert CID to integer to avoid .0 display
        try:
            cid_int = int(float(key))
            newarr.append(str(cid_int))  # Add the compound ID as integer string
        except (ValueError, TypeError):
            newarr.append(str(key))  # Fallback to original if conversion fails
            
        # === CRITICAL: Validate index bounds before accessing val elements ===
        for i in range(6 + num_containers):
            if i >= len(val):
                logger.error(f"Index {i} out of bounds for val (length {len(val)}) for CID {key}")
                newarr.append("N/A")  # Safe fallback
            else:
                newarr.append(val[i])  # Add the normalized values

        cid = int(key)
        matching_rows = df.index[df["MoleculeID"] == cid].tolist()
        if matching_rows:
            matching_row = matching_rows[0]
            for colname in checkstring:
                try:
                    if colname in df.columns:
                        newarr.append(str(df.loc[matching_row, colname]))
                    else:
                        newarr.append("N/A")  # Column doesn't exist
                except Exception as e:
                    newarr.append("N/A")  # Error accessing value
        else:
            # No matching row found, fill with N/A
            for colname in checkstring:
                newarr.append("N/A")

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
    - finalresults: The formatted final results.
    - cidarr: A list of CIDs for PDF generation.
    """
    try:
        logging.debug(f"formatResults called with:")
        logging.debug(f"finalresults type: {type(finalresults)}")
        logging.debug(f"finalresults length: {len(finalresults) if hasattr(finalresults, '__len__') else 'N/A'}")
        logging.debug(f"qcid: {qcid}")
        logging.debug(f"qcid type: {type(qcid)}")

        if finalresults and len(finalresults) > 0:
            logging.debug(f"First result: {finalresults[0]}")
            logging.debug(f"First result type: {type(finalresults[0])}")
            if hasattr(finalresults[0], '__len__'):
                logging.debug(f"First result length: {len(finalresults[0])}")

        # Handle case where qcid is None or empty
        if qcid is None or qcid == "":
            logging.debug(f"No qcid provided, returning empty results")
            return [], []
        
        # Note: qcid == -1 is valid (query not in PubChem but has OPERA predictions)

        # Limit to top 20 results for PDF generation, plus query
        if len(finalresults) > 21:
            finalresults = finalresults[:21]
            logging.debug(f"After slicing to 21: {len(finalresults)} results")

        cidarr = [str(qcid)]
        for i in range(len(finalresults)):
            logging.debug(f"Processing result {i}: {finalresults[i]}")
            # Ensure the result is not empty
            if finalresults[i] and len(finalresults[i]) > 0:
                current_cid = str(finalresults[i][0])
                logging.debug(f"Current CID: '{current_cid}', str(qcid): '{str(qcid)}'")
                
                # Add to cidarr if it's not the query and not already present
                if current_cid != str(qcid) and current_cid not in cidarr:
                    cidarr.append(current_cid)
                    logging.debug(f"Added CID to array: {current_cid}")
                else:
                    logging.debug(f"CID already in array or matches query: {current_cid}")

                # Format the first column to be just the CID for display
                if "Query" in str(finalresults[i][0]):
                    finalresults[i][0] = qcid
                    logging.debug(f"Marked query result: {finalresults[i][0]}, stored qcid: {qcid}")
                else:
                    finalresults[i][0] = finalresults[i][0]
                
                logging.debug(f"Formatted result {i}: {finalresults[i][0]}")
            else:
                logging.debug(f"Skipping empty result at index {i}")

        logging.debug(f"Final cidarr: {cidarr}")
        logging.debug(f"Final finalresults length: {len(finalresults)}")
        return finalresults, cidarr
    except Exception as e:
        logging.debug(f"Error in formatResults: {e}")
        import traceback
        logging.debug(f"Traceback: {traceback.format_exc()}")
        return finalresults, [str(qcid)]
