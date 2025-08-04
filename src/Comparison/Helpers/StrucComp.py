# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pandas as pd
from utils.progress_logger import get_progress_logger

def extraStrucComp(comp_dict, job_id="default"):
    """
    Performs additional structural comparisons and updates the comparison dictionary.
    Uses static OPERA column names.
    """
    logger = get_progress_logger(job_id)
    logger.info("Starting extra structural comparisons.")

    try:
        # Read the CSV file
        df = pd.read_csv("src/Comparison/LocalIO/Thermout.csv")
        logger.info(f"Loaded Thermout.csv with shape {df.shape}")
        
        # Use static OPERA column names
        molecule_id_column = 'MoleculeID'
        molecular_weight_column = 'MolWeight'
        
        logger.info(f"Available columns: {df.columns.tolist()}")
        
        if molecule_id_column not in df.columns:
            logger.error(f"Required column {molecule_id_column} not found in DataFrame")
            return comp_dict
        
        cids = df[molecule_id_column].values
        finalarr = [1] * len(cids)
        qcid = cids[0]
        comp_dict[str(qcid)][2] = 1
        cids = cids[1:]

        # Perform molecular weight comparison
        if molecular_weight_column in df.columns:
            logger.debug(f"Processing molecular weight using column: {molecular_weight_column}")
            vals = df[molecular_weight_column].values
            comp = vals[0]
            if float(comp) == 0:
                comp = 0.0001
            vals = vals[1:]
            for i in range(len(cids)):
                key = str(cids[i])
                if key not in comp_dict:
                    logger.warning(f"CID {key} not found in comp_dict during molecular weight update")
                    continue
                if len(comp_dict[key]) <= 2:
                    logger.warning(f"comp_dict[{key}] has insufficient elements for molecular weight: {comp_dict[key]}")
                    continue
                if i >= len(vals):
                    logger.error(f"Index {i} out of bounds for molecular weight vals (length {len(vals)})")
                    continue
                    
                finalval = 1 / (1 + abs((float(vals[i]) - comp) / (comp)))
                comp_dict[key][2] = finalval
        else:
            logger.warning("No molecular weight column found, skipping molecular weight comparison")

        # Define weights and structural properties using static OPERA column names
        weights = [0.1, 0.1, 0.1, 0.1]
        comparison_properties = [
            'nbRing',
            'nbLipinskiFailures', 
            'TopoPolSurfAir',
            'nbC'
        ]

        property_names = ['rings', 'lipinski_failures', 'topological_polar_surface_area', 'carbons']

        logger.info(f"Using structural comparison columns: {comparison_properties}")

        for j in range(len(weights)):
            column_name = comparison_properties[j]
            prop_name = property_names[j]
            
            if column_name in df.columns:
                logger.debug(f"Processing {prop_name} using column: {column_name}")
                vals = df[column_name].values
                comp = vals[0]
                vals = vals[1:]
                for i in range(len(cids)):
                    if float(comp) == 0:
                        finalarr[i] *= (1 / (1 + abs((float(vals[i]) / 0.001)))) ** weights[j]
                    else:
                        finalarr[i] *= (
                            1 / (1 + abs((float(vals[i]) - comp) / comp))
                        ) ** weights[j]
                for i in range(len(cids)):
                    key = str(cids[i])
                    if key not in comp_dict:
                        logger.warning(f"CID {key} not found in comp_dict during structural property update")
                        continue
                    if len(comp_dict[key]) <= 1:
                        logger.warning(f"comp_dict[{key}] has insufficient elements for structural update: {comp_dict[key]}")
                        continue
                    if i >= len(finalarr):
                        logger.error(f"Index {i} out of bounds for finalarr (length {len(finalarr)})")
                        continue
                        
                    comp_dict[key][1] *= finalarr[i]
            else:
                logger.warning(f"No column found for {prop_name}, skipping this comparison")

        logger.info("Extra structural comparisons completed successfully.")
        return comp_dict
        
    except Exception as e:
        logger.error(f"Error during extra structural comparisons: {e}")
        raise
