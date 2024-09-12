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

# Exception Count is a global variable that counts how many CIDs faced error and were ignored
exceptioncount = 0


# Function takes a smile and returns fingerprint, increments exception count if failing
def comp_fp(smile):
    global exceptioncount
    molec = Chem.MolFromSmiles(smile)
    if molec is not None:
        return AllChem.GetMorganFingerprintAsBitVect(
            molec, radius=2, nBits=2048, useFeatures=True
        ).ToBitString()
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
            print("chunk done")

    # #nrows underneath set to 100000, parameter should be removed when actually using, will take very long
    # df = pd.read_csv(path, nrows = 1000, delimiter = '\t', header = None, names = ['cid', 'smiles'])

    # # change smiles column to be fingerprints and save to csv
    # df["smiles"] =  df["smiles"].apply(comp_fp)
    # df.rename(columns={'smiles': 'fps'}, inplace = True)
    # df.to_csv(path2, index = False)

    print("Exception Count: " + str(exceptioncount))
    print("Elapsed time: " + str(time.time() - current_time))
