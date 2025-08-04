# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pandas as pd
from utils.progress_logger import get_progress_logger

def get_column_safely(dataframe, row_index, column_name, default_value=0):
    """
    Safely get a value from a DataFrame column with fallback.
    """
    if column_name not in dataframe.columns:
        return default_value
    
    try:
        value = dataframe.loc[row_index, column_name]
        if value != "N/A" and value is not None and str(value).strip() != "":
            return float(value)
        else:
            return default_value
    except (ValueError, TypeError, KeyError):
        return default_value


def toxicComparison(comp_dict, job_id="default"):
    """
    Performs toxicity property comparisons and updates the comparison dictionary.
    Uses static OPERA column names.
    """
    logger = get_progress_logger(job_id)
    logger.info("Starting toxicity property comparisons.")

    try:
        df = pd.read_csv("src/Comparison/LocalIO/Thermout.csv")
        logger.info(f"Loaded Thermout.csv with shape {df.shape}")
        
        # Use static OPERA column names
        molecule_id_column = 'MoleculeID'
        
        logger.info(f"Available columns: {df.columns.tolist()}")
        
        if molecule_id_column not in df.columns:
            logger.error(f"Required column {molecule_id_column} not found in DataFrame")
            return comp_dict
        
        for cid, vals in comp_dict.items():
            cid_int = int(cid)
            
            # Validate comp_dict structure before modification
            if not isinstance(vals, (list, tuple)):
                logger.warning(f"comp_dict[{cid}] is not a list/tuple: {type(vals)}")
                continue
            if len(vals) <= 4:
                logger.warning(f"comp_dict[{cid}] has insufficient elements for toxicity score: {vals}")
                continue
            
            # Find the matching row
            matching_rows = df.index[df[molecule_id_column] == cid_int].tolist()
            if not matching_rows:
                logger.warning(f"No matching row found for CID {cid}")
                continue
                
            matching_row = matching_rows[0]
            
            # Use static OPERA column names
            BCF = get_column_safely(df, matching_row, 'LogBCF_pred', default_value=1.0)
            EPA_raw = get_column_safely(df, matching_row, 'CATMoS_EPA_pred', default_value=0.0)
            LD50 = get_column_safely(df, matching_row, 'CATMoS_LD50_pred', default_value=1000.0)
            
            EPA = 1 + EPA_raw
            LD50_scaled = LD50 / 1000
            totalval = BCF * EPA * LD50_scaled / 2
            
            vals[4] = totalval
            
            logger.debug(f"CID {cid}: BCF={BCF}, EPA={EPA}, LD50={LD50_scaled}, totalval={totalval}")
            
        logger.info("Toxicity property comparisons completed successfully.")
        return comp_dict
        
    except Exception as e:
        logger.error(f"Error during toxicity property comparisons: {e}")
        raise
