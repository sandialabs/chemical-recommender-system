# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import os, sys, copy, logging
from rdkit import Chem
from rdkit.Chem import AllChem, RDConfig

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # type: ignore

from Comparison.Helpers.ThermComp import thermalComparison as tc
from Comparison.Helpers.ToxicComp import toxicComparison as toxc
from Comparison.Helpers.StrucComp import extraStrucComp as sc

from Comparison.utils.gen import parseQuery, fillDict

from Comparison.Runners.MilvusRunner import runMilvus as fp
from Comparison.Runners.OperaRunner import runOpera as rp
from Comparison.Runners.ModelRunner import runModels as rm


def setup_logging():
    """
    Set up logging configuration to refresh the log file each time this function is called.
    """
    # Clear any existing loggers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    # Configure the main log file
    comparison_handler = logging.FileHandler("logs/comparison.log", mode="w")
    comparison_handler.setLevel(logging.DEBUG)
    comparison_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    comparison_handler.setFormatter(comparison_formatter)

    # Configure the root logger to also write to a separate file
    root_handler = logging.FileHandler("logs/comparison-root.log", mode="w")
    root_handler.setLevel(logging.INFO)
    root_formatter = logging.Formatter("%(message)s")
    root_handler.setFormatter(root_formatter)

    # Get the logger for this module and the root logger
    logger = logging.getLogger(__name__)
    root_logger = logging.getLogger()

    # Add handlers to the loggers
    logger.addHandler(comparison_handler)
    root_logger.addHandler(root_handler)

    # Set specific loggers to higher levels to suppress their debug logs
    logging.getLogger("PIL.PngImagePlugin").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
    logging.getLogger("pubchempy").setLevel(logging.WARNING)


def comparisonFunction(
    queryinput,
    finnum,
    tarray,
    incele,
    smarts,
    smarts_num,
    tries=1,
    containers=None,
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

    Returns:
    - A list containing the comparison dictionary, query models, and a flag indicating if substructure searching failed.
    """
    # Save the existing loggers
    existing_handlers = logging.root.handlers[:]

    # Set up new logging configuration
    setup_logging()
    logger = logging.getLogger(__name__)
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
                query = AllChem.GetMorganFingerprintAsBitVect(
                    Chem.MolFromSmiles(query_smi),
                    radius=2,
                    nBits=2048,
                    useFeatures=True,
                )
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
        heap = fp(query, heapnum)

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

        # Run Opera Comparison from command line
        rp(comp_dict, querycid, tarray, query_smi, smiles_dict)

        logger.info("============ Analyzing OPERA Results ============")
        # Calculate dict values for each of the Opera metrics
        comp_dict = sc(comp_dict)
        comp_dict = tc(comp_dict, tarray)
        comp_dict = toxc(comp_dict)

        # For each item in the dict, compute the rdkit SA score and add this to the consideration
        logger.info("Computing SA Scores")
        for key, val in comp_dict.items():
            if key == "-1":
                smiles = query_smi
            else:
                smiles = smiles_dict[key]

            mol = Chem.MolFromSmiles(smiles)
            score = 10 - sascorer.calculateScore(mol)
            val[0] = val[0] * score
            val[5] = score

        ###########################################################################################################################################################################################################

        # Check if CLI users specified extra models to be used, add in now to comp
        first = True
        query_models = []
        if containers:
            for image in containers:
                if first:
                    logger.info("Running Added Models...")
                    first = False
                out = rm(comp_dict, image, querysmiles, smiles_dict)
                comp_dict = out[0]
                query_models.append(out[1])

        ###########################################################################################################################################################################################################

        logger.info("Comparison function completed successfully.")
        return [comp_dict, query_models, subfailed]
    except Exception as e:
        logger.error(f"Error during comparison function: {e}")
        raise
    finally:
        # Revert to the previous logging configuration
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        for handler in existing_handlers:
            logging.root.addHandler(handler)


# ComparisonFunction("1923", 10, [True,True,False,False,False], False, 'on', None, None, containers=['example'])

# ComparisonFunction("C2(CC1(=CC=C(C=C1)N))(=CC=C(C=C2)N)", 10, [True,True,False,False,False], False, 'on', '[NH2]', 2) # <- Test params

# ComparisonFunction("CCCCCCCCC1C(C=CC(C1CCCCCCCC(=O)O)CCCCCC)CCCCCCCC(=O)O", 1, [True,True,True,True,True], False, True, ' [#6R2][#8R1;r3][#6R2]', 2)

# ComparisonFunction("CCCCCCCCC2C(CCCCCCCC(=O)OCC1CO1)C=CC(CCCCCC)C2CCCCCCCC(=O)OCC3CO3", 10,  [True,True,True,True,True], False, True, None, None) <- no cid
