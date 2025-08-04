# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import os
import sys
import pandas as pd
import pubchempy as pcp
from flask import current_app

# Add the OPERA directory to the Python path
opera_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'OPERA')
sys.path.append(opera_path)

from opera_cli import OPERACLI

# Import centralized logging
from utils.progress_logger import get_progress_logger

# Get a logger for this module (will be replaced with job-specific logger)
logger = None

def opera_predictions_to_dataframe(predictions_data, cid_list):
    """
    Convert OPERA prediction results to a pandas DataFrame.
    
    Args:
        predictions_data: List of dictionaries from OPERA predictions
        cid_list: List of molecule IDs corresponding to the predictions
        
    Returns:
        DataFrame of OPERA predictions
    """
    if not predictions_data:
        return pd.DataFrame()
    
    # Create DataFrame from predictions
    df = pd.DataFrame(predictions_data)
    
    # Ensure MoleculeID column matches the input CID list
    if 'MoleculeID' in df.columns:
        # Replace MoleculeID with actual CIDs if they don't match
        if len(cid_list) == len(df):
            df['MoleculeID'] = cid_list
    else:
        # Add MoleculeID if missing
        df['MoleculeID'] = cid_list[:len(df)]
    
    return df

def smiles_to_opera_properties(smiles_list, cid_list, job_id="default", requested_properties=None):
    """
    Generate molecular properties using OPERA.
    
    Args:
        smiles_list: List of SMILES strings
        cid_list: List of molecule IDs corresponding to the SMILES
        job_id: Unique identifier for this job (for logging purposes)
        requested_properties: List of property keys to compute
    
    Returns:
        bool: True if OPERA succeeded, False if OPERA failed
    """
    # Get the app state (singleton or passed in)
    try:
        opera_prop_cache = current_app.state.opera_prop_cache
    except:
        opera_prop_cache = None

    # Build cache key from all relevant inputs
    cache_key = (
        tuple(smiles_list),
        tuple(cid_list),
        tuple(requested_properties) if requested_properties else None
    )
    cache_path = "src/Comparison/LocalIO/Thermout.csv"

    # Check cache before running OPERA
    if opera_prop_cache and cache_key in opera_prop_cache:
        logger = get_progress_logger(job_id)
        logger.info("Returning cached OPERA property predictions.")
        opera_prop_cache.move_to_end(cache_key)
        cached_df = opera_prop_cache[cache_key]
        cached_df.to_csv(cache_path, index=False)
        return True
    logger = get_progress_logger(job_id)
    logger.progress("Processing", f"Computing OPERA properties for {len(smiles_list)} compounds")
    
    try:
        # Initialize OPERA CLI
        opera_cli = OPERACLI()
        # Send progress update for larger datasets
        if len(smiles_list) > 5:
            # Read predictions from Thermout.csv and cache them
            try:
                df = pd.read_csv(cache_path)
                if opera_prop_cache is not None:
                    opera_prop_cache[cache_key] = df
                    if len(opera_prop_cache) > 10:
                        opera_prop_cache.popitem(last=False)
            except Exception as cache_e:
                logger.warning(f"Could not cache OPERA predictions: {cache_e}")
            logger.progress("Processing", f"Running OPERA predictions for {len(smiles_list)} compounds...")
        # Run OPERA predictions with requested properties
        opera_success = opera_cli.run_opera(smiles_list, cid_list, requested_properties)
        if opera_success:
            # OPERA succeeded - predictions are written to Thermout.csv
            logger.info("OPERA CLI execution completed successfully")
            return True
        else:
            logger.warning("OPERA CLI execution failed - attempting fallback per-compound OPERA runs")
            # Fallback: try each compound individually
            successful_rows = []
            failed_cids = []
            for smiles, cid in zip(smiles_list, cid_list):
                try:
                    logger.progress("Processing", f"Running OPERA predictions for CID {cid}...")
                    single_success = opera_cli.run_opera([smiles], [cid], requested_properties)
                    logger.info("Processing", f"OPERA predictions for CID {cid} completed.")
                    logger.info(f"OPERA predictions for CID {cid} {'succeeded' if single_success else 'failed'}")
                    logger.info(f"{single_success}")
                    if single_success:
                        # Read the result from Thermout.csv and append to successful_rows
                        df = pd.read_csv(cache_path)
                        df['MoleculeID'] = df['MoleculeID'].astype(str)
                        cid_str = str(cid)
                        row = df[df['MoleculeID'] == cid_str]
                        if row.isnull().values.any():
                            logger.warning(f"Retrying OPERA for CID {cid} due to missing predictions.")
                            retry = opera_cli.run_opera([smiles], [cid], requested_properties)
                            if retry:
                                df = pd.read_csv(cache_path)
                                df['MoleculeID'] = df['MoleculeID'].astype(str)
                                row = df[df['MoleculeID'] == cid_str]
                        if not row.empty and not row.isnull().values.any():
                            successful_rows.append(row.iloc[0].to_dict())
                        else:
                            failed_cids.append(cid)
                    else:
                        failed_cids.append(cid)
                except Exception as indiv_e:
                    logger.warning(f"OPERA failed for CID {cid}: {indiv_e}")
                    failed_cids.append(cid)
            if successful_rows:
                # Write all successful rows to Thermout.csv
                pd.DataFrame(successful_rows).to_csv(cache_path, index=False)
                # Cache the fallback results once after all runs
                try:
                    df = pd.read_csv(cache_path)
                    if opera_prop_cache is not None:
                        opera_prop_cache[cache_key] = df
                        if len(opera_prop_cache) > 10:
                            opera_prop_cache.popitem(last=False)
                except Exception as cache_e:
                    logger.warning(f"Could not cache fallback OPERA predictions: {cache_e}")
                logger.info(f"OPERA fallback succeeded for {len(successful_rows)} compounds; failed for {failed_cids}.")
                if failed_cids:
                    logger.warning(f"CIDs with failed OPERA predictions: {failed_cids}")
                # The caller can remove failed_cids from comp_dict/smiles_dict as needed
                return True
            else:
                logger.error("All fallback OPERA runs failed. No valid predictions.")
                return False
    except Exception as e:
        logger.error(f"OPERA FAILED: {e}")
        logger.progress("Warning", f"OPERA property calculation failed: {str(e)}")
        return False
    
def runProperty(comp_dict, querycid, tarray, query_smi, smiles_dict, job_id="default"):
    """
    Compute properties using OPERA and write to Thermout.csv.
    
    Args:
        comp_dict: Dictionary of candidate compounds
        querycid: Query compound ID
        tarray: Array of boolean values indicating which properties to calculate
                [MP, BP, LogP, HLaw, VP] - plus always include BCF and CATMoS
        query_smi: Query SMILES string
        smiles_dict: Dictionary mapping CIDs to SMILES
        job_id: Unique identifier for this job
        
    Returns:
        bool: True if OPERA succeeded, False if OPERA failed
    """
    logger = get_progress_logger(job_id)
    logger.info("============ Entering OPERA Property Models ============")
    
    # Create directory if it doesn't exist
    os.makedirs("src/Comparison/LocalIO", exist_ok=True)
    
    # Map tarray to property names based on the system's expectations
    # tarray[0] = MP, tarray[1] = BP, tarray[2] = LogP, tarray[3] = HLaw, tarray[4] = VP
    requested_properties = []
    if tarray[0]: requested_properties.append('-MP')
    if tarray[1]: requested_properties.append('-BP')
    if tarray[2]: requested_properties.append('-LogP')
    if tarray[3]: requested_properties.append('-HL')
    if tarray[4]: requested_properties.append('-VP')

    requested_properties.extend(['-BCF', '-CATMoS', '-StrP'])  # Always include BCF and CATMoS

    logger.info(f"Requested OPERA properties: {requested_properties}")
    
    # Prepare SMILES and CID arrays
    if query_smi:
        querysmiles = query_smi
    else:
        try:
            querysmiles = pcp.get_properties("SMILES", querycid)[0]["SMILES"]
        except Exception as e:
            logger.error(f"Error getting SMILES for CID {querycid}: {e}")
            querysmiles = ""
    
    # Build arrays for all compounds (query + candidates)
    smilesarr = [querysmiles]
    cidarr = [querycid]
    
    for key, _ in comp_dict.items():
        cidarr.append(key)
        if key == "-1":
            smilesarr.append(query_smi)
        else:
            smilesarr.append(smiles_dict.get(key, ""))
    
    # Filter out empty SMILES
    valid_pairs = [(smiles, cid) for smiles, cid in zip(smilesarr, cidarr) if smiles and smiles.strip()]
    
    if not valid_pairs:
        logger.error("No valid SMILES found for property calculation")
        return False
    
    valid_smiles, valid_cids = zip(*valid_pairs)
    
    # Run OPERA - it either succeeds (and creates Thermout.csv) or fails
    opera_success = smiles_to_opera_properties(list(valid_smiles), list(valid_cids), job_id, requested_properties)
    
    if opera_success:
        logger.info("OPERA property predictions completed successfully")
        return True
    else:
        logger.warning("OPERA property calculation failed")
        return False
