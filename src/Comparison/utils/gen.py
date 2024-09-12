# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import heapq
import re
import time
import ast
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp

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


def get_compound_properties(cid):
    """
    Get necessary properties for the compound from its CID.

    Inputs:
    - cid: The chemical identifier (CID).

    Returns:
    - smiles: The canonical SMILES string of the compound.
    - queryname: The IUPAC name or a synonym of the compound.
    """
    logger.debug(f"Fetching properties for CID: {cid}")
    properties = pcp.get_properties(["CanonicalSMILES", "IUPACName"], cid)[0]
    smiles = properties.get("CanonicalSMILES", -1)
    queryname = properties.get("IUPACName", None)
    if queryname is None:
        synonyms = pcp.get_synonyms(cid)[0].get("Synonym", [])
        queryname = synonyms[0] if synonyms else f"CID: {cid}"
    logger.debug(f"Properties for CID {cid}: SMILES={smiles}, Name={queryname}")
    return smiles, queryname


def parseQuery(queryinput):
    """
    Parse the query input to determine the query fingerprint, CID, name, and SMILES.

    Inputs:
    - queryinput: The input query, which can be a CID, name, or SMILES string.

    Returns:
    - A list containing the query fingerprint, query CID, query name, and query SMILES.
    """
    logger.debug(f"Parsing query input: {queryinput}")
    # Start with bad values for parity at the end
    query, querycid, queryname, smiles = -1, -1, -1, -1

    # Try searching by CID first, assume so if input is all numbers
    if isinstance(queryinput, int) or queryinput.isdigit():
        querycid = int(queryinput)
        smiles, queryname = get_compound_properties(querycid)
        if smiles != -1:
            query = AllChem.GetMorganFingerprintAsBitVect(
                Chem.MolFromSmiles(smiles), radius=2, nBits=2048, useFeatures=True
            )
    # Next, attempt to search as from name
    if query == -1:
        try:
            querycid = pcp.get_cids(queryinput, "name", list_return="flat")[0]
            smiles, queryname = get_compound_properties(querycid)
            if smiles != -1:
                query = AllChem.GetMorganFingerprintAsBitVect(
                    Chem.MolFromSmiles(smiles), radius=2, nBits=2048, useFeatures=True
                )
        except:
            pass

    # Next, attempt to search as from smiles
    if query == -1:
        try:
            querycid = pcp.get_compounds(queryinput, "smiles")[0].cid
            smiles, queryname = get_compound_properties(querycid)
            if smiles != -1:
                query = AllChem.GetMorganFingerprintAsBitVect(
                    Chem.MolFromSmiles(smiles), radius=2, nBits=2048, useFeatures=True
                )
        except:
            pass

    logger.debug(
        f"Parsed query result: {query}, CID: {querycid}, Name: {queryname}, SMILES: {smiles}"
    )
    return [query, querycid, queryname, smiles]


def filterCIDs(found, queryname, incEle, smarts_mol, smarts_num, querysmiles):
    """
    Filter the list of CIDs based on user parameters and return a list of passing SMILES strings.

    Inputs:
    - found: A list of tuples containing the CIDs and their similarity scores.
    - queryname: The name of the query compound.
    - incEle: A list of included elements or True to include all elements.
    - smarts_mol: The SMARTS pattern for substructure searching.
    - smarts_num: The number of substructure matches required.
    - querysmiles: The SMILES string of the query compound.

    Returns:
    - c_arr: A list of SMILES strings for the filtered CIDs.
    """
    logger.debug("Filtering CIDs based on user parameters.")
    allowed = ["H", "C", "N", "O", "F", "P", "S", "Cl", "Se", "Br", "I"]

    if incEle != True:
        for element in incEle:
            allowed.append(element)
    # c_arr is the return parameter. For the list of CIDs inputted, if they pass the filter we return its SMILES and otherwise None. c_arr is ordered the same as cids.
    c_arr = []
    cids = []
    for item in found:
        cids.append(item[1])

    properties = pcp.get_properties(
        ["CanonicalSMILES", "IsomericSMILES", "IUPACName", "MolecularFormula"], cids
    )

    for i in range(len(cids)):
        property = properties[i]
        cid = cids[i]
        # Extract the properties
        canonical_smiles = property.get("CanonicalSMILES", "")
        isomeric_smiles = property.get("IsomericSMILES", "")
        iupac_name = property.get("IUPACName", "")
        elements = re.findall(r"[A-Z][a-z]*", property.get("MolecularFormula", ""))

        # If the name is not defined by IUPAC, search for existing synonyms
        name = iupac_name
        if name is None:
            try:
                synonyms = pcp.get_synonyms(cid)[0]["Synonym"]
                if synonyms:
                    name = synonyms[0]
                else:
                    name = "CID: " + str(cid)
            except:
                name = "CID: " + str(cid)

        name1 = queryname + ";"
        name2 = ";" + queryname
        name3 = queryname + ","
        name4 = "," + queryname

        # Filter out results that contain the original query value as a part of the molecule
        if (
            name is None
            or queryname is None
            or name1 in name
            or name2 in name
            or name3 in name
            or name4 in name
            or ";" in name
            or name is queryname
        ):
            c_arr.append(None)
            continue

        # Filter out non-bonded molecules
        if "." in str(canonical_smiles):
            c_arr.append(None)
            continue

        # Check if molecule only has allowed elements
        if incEle != True:
            found = False
            elements = elements
            for i in elements:
                if i not in allowed:
                    c_arr.append(None)
                    found = True
                    continue
            if found is True:
                continue

        # Perform substructure searching if specified by user
        if smarts_mol:
            if not smarts_num:
                if not (
                    Chem.MolFromSmiles(isomeric_smiles).GetSubstructMatches(smarts_mol)
                ):
                    c_arr.append(None)
                    continue
            else:
                matches = Chem.MolFromSmiles(isomeric_smiles).GetSubstructMatches(
                    smarts_mol
                )
                if len(matches) != int(smarts_num):
                    c_arr.append(None)
                    continue

        # If passing all the above tests, add the molecule SMILES to c_arr for return

        if canonical_smiles == querysmiles:
            c_arr.append(None)
            continue
        c_arr.append(canonical_smiles)
    logger.debug(f"Filtered CIDs: {c_arr}")
    return c_arr


def fillDict(
    comp_dict, heap, finsize, queryname, incEle, smarts_mol, smarts_num, querysmiles
):
    """
    Filter candidates to create a shortlist and fill the comparison dictionary.

    Inputs:
    - comp_dict: A dictionary where keys are chemical identifiers (CIDs) and values are lists containing comparison metrics.
    - heap: A heap containing the CIDs and their similarity scores.
    - finsize: The final size of the comparison dictionary.
    - queryname: The name of the query compound.
    - incEle: A list of included elements or True to include all elements.
    - smarts_mol: The SMARTS pattern for substructure searching.
    - smarts_num: The number of substructure matches required.
    - querysmiles: The SMILES string of the query compound.

    Returns:
    - smiles_dict: A dictionary mapping CIDs to SMILES strings for the filtered candidates.
    """
    logger.info("============= Filtering candidates for shortlist ============")
    filter_start = time.time()
    firstpass = True
    smiles_dict = {}
    while len(comp_dict) < finsize:
        # Takes the heap and changes the highest values of finsize into dict format
        # Checks it through the filter cid function and if it works, the cid is added to the dict
        # Dict is of the form {key: [similarity, similarity, ...]}. 0th parameter is updated as through program as aggregate similarity while 1st is the constant similarity for future reference.
        if firstpass:
            found = heapq.nlargest(finsize, heap)
            firstpass = False
        else:
            found = heapq.nlargest(int(finsize / 2), heap)
        if len(found) == 0:
            break

        vals = filterCIDs(found, queryname, incEle, smarts_mol, smarts_num, querysmiles)
        for i in range(len(found)):
            item = found[i]
            heap.remove(item)
            if vals[i] is not None:
                if len(comp_dict) < finsize:
                    smiles_dict[item[1]] = vals[i]
                    comp_dict[item[1]] = [item[0], item[0], None, None, None, None]
    filter_end = time.time()
    logger.info(
        f"Total time taken for filtering: {filter_end - filter_start:.2f} seconds"
    )

    return smiles_dict


def ensure_quoted_elements(element_list):
    """
    Ensure that all elements in the list are properly quoted strings.
    """
    if isinstance(element_list, str):
        element_list = element_list.strip("[]").split(",")
    return [
        (
            f'"{elem.strip()}"'
            if not (elem.startswith('"') and elem.endswith('"'))
            else elem.strip()
        )
        for elem in element_list
    ]


def parseBatchInput(batch_text, from_command=False, containers=[]):
    """
    Parse the batch input text and create elements for each query.

    Inputs:
    - batch_text: The input text containing multiple queries.
    - from_command: A flag indicating whether the input is from a command line text file.

    Returns:
    - queries: A list of dictionaries containing parsed query parameters.
    """
    logger.debug("Parsing batch input text.")
    queries = []

    # Slightly different parsing technique if submitted via command line text file or through GUI
    if from_command:
        lines = batch_text.split("\n")
    else:
        lines = batch_text.split("\r\n")

    # Evaluate the lines, pull in data into each value as needed
    # Input format: query, num candidates, [Melting Point, Boiling Point, logP, Henry's Law, Vapor Pressure], Include All Elements, Substructure Smarts, Substructure Hits, Weightages (optional)
    # Split by commas into array and parse each array element accordingly
    for line in lines:
        queryinput = line.split(", ")
        query = queryinput[0]
        if query == "":
            continue

        try:
            finnum = int(queryinput[1])
            logger.debug(f"finnum: {finnum}")

            Tstring = queryinput[2][1:-1].split(",")
            tarray = [i == "True" for i in Tstring]
            logger.debug(f"tarray: {tarray}")

            incEle = queryinput[3]
            logger.debug(f"incEle (initial): {incEle}")

            if incEle == "True":
                incEle = True
            else:
                if queryinput[4] != "None":
                    incEle = ensure_quoted_elements(queryinput[4])
                    incEle = ast.literal_eval(f"[{', '.join(incEle)}]")
                else:
                    incEle = [""]
            logger.debug(f"incEle (processed): {incEle}")

            smarts = queryinput[5]
            if smarts == "None":
                smarts = None
            logger.debug(f"smarts: {smarts}")

            smarts_num = int(queryinput[6]) if queryinput[6] != "None" else None
            logger.debug(f"smarts_num: {smarts_num}")

            # Only look for extra weightage in the input if there exist extra parameters
            # The optional weightages should be inputted in the format [Structural Similarity, Molecular Weight, Thermophysical Similarity, Predicted Toxicity, SA Scoring]
            if len(queryinput) < 8:
                weights = [1] * (5 + len(containers))
            else:
                weights = ast.literal_eval(queryinput[7])
                if len(weights) != (len(containers) + 5):
                    raise ValueError(
                        "Weightage array is not of correct length  - ensure added models are accounted for."
                    )
            logger.debug(f"weights: {weights}")

            queries.append(
                {
                    "query": query,
                    "finnum": finnum,
                    "tarray": tarray,
                    "incEle": incEle,
                    "smarts": smarts,
                    "smarts_num": smarts_num,
                    "weights": weights,
                }
            )
        except (ValueError, SyntaxError) as e:
            logger.error(f"Error parsing line: {line}. Error: {e}")
            continue

    logger.debug(f"Parsed batch queries: {queries}")
    return queries
