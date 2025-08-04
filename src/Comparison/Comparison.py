# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import os, sys, copy, logging, pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, RDConfig

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # type: ignore

from Comparison.Helpers.ThermComp import thermalComparison as tc
from Comparison.Helpers.ToxicComp import toxicComparison as toxc
from Comparison.Helpers.StrucComp import extraStrucComp as sc

from Comparison.utils.gen import parseQuery, fillDict

from Comparison.Runners.MilvusRunner import runMilvus as fp
from Comparison.Runners.OperaRunner import runProperty as rp
from Comparison.Runners.ModelRunner import runModels as rm

# Import the unified logging
from utils.progress_logger import get_progress_logger

# Get a logger for this module (will be replaced with job-specific logger in function)
logger = None

# Function to be called at the start of a comparison run
def setup_logging_for_comparison(job_id="default"):
    """
    Set up logging for a new comparison run, using the unified progress logger.
    """
    return get_progress_logger(job_id)


def comparisonFunction(
    queryinput,
    finnum,
    tarray,
    incele,
    smarts,
    smarts_num,
    tries=1,
    containers=None,
    job_id="default",
    disallow_isotopes=False,
):
    """
    Perform a comprehensive comparison of chemical compounds based on various metrics.

    Inputs:
    - queryinput: The input query, which can be a CID, name, or SMILES string.
    - finnum: The number of final results to be returned.
    - tarray: An array of boolean values indicating which thermophysical properties are to be calculated.
      The order is:
        - Melting Point (MP)
        - Boiling Point (BP)
        - Log P (logP)
        - Vapor Pressure (VP)
        - Henry's Law Constant (Hlaw)
    - incele: A flag indicating whether to include all elements.
    - smarts: The SMARTS pattern for substructure searching.
    - smarts_num: The number of substructure matches required.
    - tries: The number of attempts to find suitable candidates (default is 1).
    - containers: A list of container names for additional models (default is None).
    - job_id: Unique identifier for this job (for logging purposes, default is "default").

    Returns:
    - A list containing the comparison dictionary, query models, and a flag indicating if substructure searching failed.
    """
    # Save the existing loggers
    existing_handlers = logging.root.handlers[:]

    # Set up new logging configuration with job-specific logger
    logger = get_progress_logger(job_id)
    logger.info("============ Starting comparison function ============")

    try:
        # Set the smiles to be queried for, use helper function to determine method of input and set variables
        query, querycid, queryname, querysmiles = parseQuery(queryinput)
        logger.debug(
            f"Parsed query: {query}, CID: {querycid}, Name: {queryname}, SMILES: {querysmiles}"
        )

        # query_smi is used in cases where the input query is actually just SMILES unknown to PubChem
        query_smi = None

        # Check if the query exists in Pubchem, else we reference its CID as -1 and use its SMILES representation. Create a FP now to move forward
        if query == -1 or querycid == -1:
            if Chem.MolFromSmiles(queryinput) is not None:
                logger.warning(
                    f"Molecule not found in PubChem, searching by SMILES: {queryinput}"
                )
                query_smi = queryinput
                querysmiles = query_smi
                querycid = -1
                queryname = ""
                fpgen = AllChem.GetMorganGenerator(
                    radius=2,
                    countSimulation=False,
                    includeChirality=False,
                    useBondTypes=True,
                    includeRingMembership=True,
                    fpSize=2048 
                )
                query = fpgen.GetFingerprint(Chem.MolFromSmiles(query_smi))

            else:
                logger.error("There is an error with the input")
                sys.exit(1)
        else:
            logger.info(f"Searching for {queryname} (cid: {querycid})")

        # Initialize variables and set start time, heapnum is how many of the top candidates from fingerprinting are kept
        # Keeping a larger number may improve accuracy while increasing runtime
        # smarts and smarts mol are information on the substructure searching

        heap = []
        smarts_mol = None
        if smarts is not None:
            smarts_mol = Chem.MolFromSmarts(smarts)

        try:
            finsize = finnum * 3
        except Exception as e:
            logger.error(f"Settings are bad input: {e}")
            sys.exit(1)

        heapnum = finsize * 20
        if tries == 2:
            heapnum = heapnum * 5
        if tries == 3:
            heapnum = heapnum * 10

        # Milvus max limit
        if heapnum > 16384:
            heapnum = 16384

        # This function will use the Milvus Vector DB in other containers to search for nearby fingerprints
        heap = fp(query, heapnum, job_id)
        
        # Debug: Check heap structure for tie handling
        if heap and len(heap) > 0:
            logger.debug(f"Heap returned {len(heap)} candidates")
            logger.debug(f"Sample heap items: {heap[:3] if len(heap) >= 3 else heap}")
            if len(heap) > 1:
                logger.debug(f"Score range: {heap[0][0]:.6f} to {heap[-1][0]:.6f}")
        heap_copy = copy.deepcopy(heap)
        comp_dict = {}

        # This function does filtering and fills up comp_dict, which is the method of keeping track of comparison metrics through this file
        # This will also return a smiles_dict that keeps track of computed SMILES for candidates to reduce recalculation latency.
        smiles_dict = fillDict(
            comp_dict,
            heap,
            finsize,
            queryname,
            incele,
            smarts_mol,
            smarts_num,
            querysmiles,
            job_id,
            disallow_isotopes=disallow_isotopes,
        )

        subfailed = False
        # If no results are found, try once more without substructure filtering, recursive call of same function with larger heap size, max 3 tries
        if len(comp_dict) == 0:
            if smarts is not None:
                if tries < 3:
                    return comparisonFunction(
                        queryinput,
                        finnum,
                        tarray,
                        incele,
                        smarts,
                        smarts_num,
                        tries + 1,
                        containers,
                        job_id,
                        disallow_isotopes=disallow_isotopes,
                    )
                subfailed = True
                # After 3 recursive calls, if still too tight, then just remove substructure requirement and use the copy of the heap now.
                logger.info(
                    "Adjusting search to remove substructure requirements: Too strict"
                )
                smarts = None
                smarts_mol = None
                smiles_dict = fillDict(
                    comp_dict,
                    heap_copy,
                    finsize,
                    queryname,
                    incele,
                    smarts_mol,
                    smarts_num,
                    querysmiles,
                    job_id,
                    disallow_isotopes=disallow_isotopes,
                )

        # If still no candidates are in the dict, params are too strict, exit
        if len(comp_dict) == 0:
            logger.error("Input parameters are too strict, please adjust and try again")
            sys.exit(1)

        ###########################################################################################################################################################################################################
        # Initialize query values in dict
        if str(querycid) not in comp_dict.keys():
            comp_dict[str(querycid)] = [1, 1, None, None, None, None]
        smiles_dict[str(querycid)] = querysmiles

        # Run Property Comparison from command line
        opera_success = rp(comp_dict, querycid, tarray, query_smi, smiles_dict, job_id)
        
        # Track OPERA failure for reporting
        opera_failed = not opera_success
        if opera_failed:
            logger.warning("OPERA failed - thermal and toxicity comparisons will use neutral values")
            logger.progress("Warning", "Property predictions failed - using neutral values for thermal/toxicity analysis")

        logger.info("============ Analyzing Property Results ============")
        
        df = pd.read_csv("src/Comparison/LocalIO/Thermout.csv")
        molecules = df['MoleculeID'].astype(int).tolist()
        df['MoleculeID'] = df['MoleculeID'].astype(str)

        # Avoid changing dict size during iteration by iterating over a static list
        for cid in list(comp_dict.keys()):
            if int(cid) not in molecules:
                logger.warning(f"CID {cid} not found in predictions, removing from consideration")
                comp_dict.pop(cid, None)
                continue

            row = df[df['MoleculeID'] == str(cid)]
            if row.isnull().values.any():
                logger.debug(list(row.values))
                logger.warning(f"CID {cid} has None or NaN values in predictions, removing from consideration")
                comp_dict.pop(cid, None)

        # Calculate dict values for each of the property metrics
        comp_dict = sc(comp_dict, job_id)
        comp_dict = tc(comp_dict, tarray, job_id)
        comp_dict = toxc(comp_dict, job_id)

        # For each item in the dict, compute the rdkit SA score and add this to the consideration
        logger.info("Computing SA Scores")
        for key, val in comp_dict.items():
            # === VALIDATION: Check val structure before accessing indices ===
            if not isinstance(val, (list, tuple)):
                logger.error(f"comp_dict[{key}] is not a list/tuple: {type(val)}")
                continue
            if len(val) < 6:
                logger.error(f"comp_dict[{key}] has insufficient elements for SA score: {val}")
                continue
                
            if key == "-1":
                smiles = query_smi
            else:
                if key not in smiles_dict:
                    logger.warning(f"Key {key} not found in smiles_dict, skipping SA score")
                    continue
                smiles = smiles_dict[key]

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    logger.warning(f"Invalid SMILES for CID {key}: {smiles}")
                    continue
                score = 10 - sascorer.calculateScore(mol)
                val[0] = val[0] * score
                val[5] = score
            except Exception as e:
                logger.warning(f"Error calculating SA score for CID {key}: {e}")
                continue

        ###########################################################################################################################################################################################################

        # Check if CLI users specified extra models to be used, add in now to comp
        first = True
        query_models = []
        if containers:
            for image in containers:
                if first:
                    logger.info("Running Added Models...")
                    first = False
                out = rm(comp_dict, image, querysmiles, smiles_dict, job_id)
                comp_dict = out[0]
                query_models.append(out[1])

        ###########################################################################################################################################################################################################

        logger.info("Comparison function completed successfully.")
        return [comp_dict, query_models, subfailed, opera_failed]
    except Exception as e:
        logger.error(f"Error during comparison function: {e}")
        raise
    finally:
        # Revert to the previous logging configuration
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        for handler in existing_handlers:
            logging.root.addHandler(handler)
