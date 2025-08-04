# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pandas as pd
from utils.progress_logger import get_progress_logger

def thermalComparison(comp_dict, tarray, job_id="default"):
    """
    Performs thermal property comparisons and updates the comparison dictionary.
    Uses static OPERA column names.
    """
    logger = get_progress_logger(job_id)
    logger.info("Starting thermal property comparisons.")

    try:
        # Define which thermophysical comparisons to do based on user input
        MP = tarray[0]
        BP = tarray[1]
        logP = tarray[2]
        Hlaw = tarray[3]
        VP = tarray[4]

        # Read Thermout.csv
        df = pd.read_csv("src/Comparison/LocalIO/Thermout.csv")
        logger.info(f"Loaded Thermout.csv with shape {df.shape}")
        
        # Use static OPERA column names
        molecule_id_column = 'MoleculeID'
        
        logger.info(f"Available columns: {df.columns.tolist()}")
        
        if molecule_id_column not in df.columns:
            logger.error(f"Required column {molecule_id_column} not found in DataFrame")
            return comp_dict
        
        cids = df[molecule_id_column].values
        finalarr = [1] * len(cids)
        cids = cids[1:]  # Skip the query (first row)

        # Define endpoints with static OPERA column names
        endpoints = [
            [MP, 'MP_pred'],
            [BP, 'BP_pred'], 
            [logP, 'LogP_pred'],
            [Hlaw, 'LogHL_pred'],
            [VP, 'LogVP_pred'],
        ]
        
        logger.info(f"Using endpoints: {[(enabled, col) for enabled, col in endpoints if enabled]}")

        for item in endpoints:
            enabled, column_name = item
            if enabled and column_name is not None:
                logger.debug(f"Processing {column_name} column")
                
                if column_name not in df.columns:
                    logger.warning(f"Column {column_name} not found in DataFrame, skipping")
                    continue
                    
                vals = df[column_name].values
                comp = vals[0]  # Query value (first row)
                
                try:
                    comp = float(comp)
                except Exception:
                    logger.warning(f"Could not convert query value for {column_name}: {comp}")
                    continue
                    
                if comp == 0:
                    comp = 0.001  # Avoid division by zero
                    
                vals = vals[1:]  # Candidate values (skip query)
                
                for i in range(len(cids)):
                    try:
                        v = float(vals[i])
                        finalarr[i] *= 1 / (1 + abs((v - comp) / comp))
                    except Exception:
                        continue
            elif enabled:
                logger.warning(f"Property requested but column not found for: {item}")

        # Update comparison dictionary with thermal scores
        for i in range(len(cids)):
            cid_str = str(cids[i])
            if cid_str not in comp_dict:
                logger.warning(f"CID {cid_str} not found in comp_dict during thermal update")
                continue
            if len(comp_dict[cid_str]) <= 3:
                logger.warning(f"comp_dict[{cid_str}] has insufficient elements for thermal score: {comp_dict[cid_str]}")
                continue
            if i >= len(finalarr):
                logger.error(f"Index {i} out of bounds for finalarr (length {len(finalarr)})")
                continue
                
            comp_dict[cid_str][3] = finalarr[i]

        logger.info("Thermal property comparisons completed successfully.")
        return comp_dict
        
    except Exception as e:
        logger.error(f"Error during thermal property comparisons: {e}")
        raise
