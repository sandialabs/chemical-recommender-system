# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

# Starts with download from PubChem, taking in a txt of CID-Smiles
# Outputs a csv of CID-Fingerprints that can be used for similarity search in later steps

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
import time
import pandas as pd
import os
import logging

# Exception Count is a global variable that counts how many CIDs faced error and were ignored
exceptioncount = 0


# Function takes a smile and returns fingerprint, increments exception count if failing
def comp_fp(smile):
    global exceptioncount
    molec = Chem.MolFromSmiles(smile)
    fpgen = AllChem.GetMorganGenerator(
        radius=2,
        countSimulation=False,
        includeChirality=False,
        useBondTypes=True,
    includeRingMembership=True,
    fpSize=2048 
    )
    
    if molec is not None:
        return fpgen.GetFingerprint(molec).ToBitString()
    else:
        exceptioncount += 1
        return None


if __name__ == "__main__":
    RDLogger.DisableLog("rdApp.*")
    # stops error log from printing, which slows down performance

    # paths to local directories for downloaded PubChem database and output file
    home = os.path.abspath(__file__)
    fpath = r"CID-SMILES.txt"
    fpath2 = r"CID-Fingerprints2.csv"

    # path = os.path.join(os.path.dirname(home), fpath)
    # path2 = ('/pscratch/panair/mol_fing/'+ fpath2)

    path = r"C:\Users\panair\Documents\PubChem Fingerprint Data\CID-SMILES_Folder\CID-SMILES.txt"
    path2 = r"C:\Users\panair\Documents\PubChem Fingerprint Data\CID-SMILES_Folder\CID-Fingerprints2.csv"

    current_time = time.time()

    #nrows below set to 100000, parameter should be removed when actually using, will take very long
    with open(path2, "w", newline="") as output_file:
        for chunk in pd.read_csv(
            path,
            nrows=10000,
            chunksize=1000,
            delimiter="\t",
            header=None,
            names=["cid", "fps"],
        ):
            chunk["fps"] = chunk["fps"].apply(comp_fp)
            chunk.to_csv(
                output_file, mode="a", index=False, header=not output_file.tell()
            )
            logging.info("chunk done")

    logging.info("Exception Count: " + str(exceptioncount))
    logging.info("Elapsed time: " + str(time.time() - current_time))
