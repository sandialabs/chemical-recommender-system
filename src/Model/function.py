# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

from rdkit import Chem
from rdkit.Chem import Descriptors


def process_smiles(smiles: str) -> float:
    """
    Takes a SMILES string as input and returns the estimated molecular weight.

    Parameters:
    smiles (str): A SMILES string representing a chemical molecule.

    Returns:
    float: The estimated molecular weight of the molecule.
    """
    try:
        # Convert the SMILES string to a molecule object
        molecule = Chem.MolFromSmiles(smiles)

        if molecule is None:
            raise ValueError("Invalid SMILES string")

        # Calculate the molecular weight
        molecular_weight = Descriptors.MolWt(molecule)

        return molecular_weight

    except Exception as e:
        print(f"An error occurred: {e}")
        return float("nan")


# Example usage
if __name__ == "__main__":
    smiles_example = "CCO"  # Ethanol
    result = process_smiles(smiles_example)
    print(f"The estimated molecular weight for {smiles_example} is {result:.2f}")
