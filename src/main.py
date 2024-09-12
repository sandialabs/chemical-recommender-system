# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import argparse, sys
from argparse import RawDescriptionHelpFormatter
from datetime import datetime


sys.path.append("Comparison")
sys.path.append("App")

from Comparison.Controller import BatchRun as br
from App.App import runFlask as rf


class CRS:
    def __init__(self):
        """
        Initialize the CRS class, set up the argument parser, and define the command-line arguments.
        """
        self.parser = argparse.ArgumentParser(
            description="""CRS - Command-line tool for managing and using the CRS

Input File Format:

Input search parameters for the run must be written into a separate file. You can use text editors like nano or vim. The format of the input should be as follows:

    query, final_number, thermo_array, include_all_elements, include_specific_elements, substructure_search, number_substructure_search

Parameter Descriptions:

- query: The query to be searched for. This can be a PubChem CID, IUPAC Name, or SMILES.
- final_number: The number of resultant candidates to be outputted in the final report.
- thermo_array: An array of 5 boolean values (True or False) to decide if the user wants to include the following thermophysical properties (in this order, no spaces between the commas):
  - Melting Point
  - Boiling Point
  - Log P
  - Vapor Pressure
  - Henry's Law Constant
- include_all_elements: CRS by default only searches these elements: H, C, N, O, F, P, S, Cl, Se, Br, I. Setting this parameter to True will include all elements instead. (Must be True or False)
- include_specific_elements: Add specific elements to search in addition to the default ones. Should be in a comma-separated format with no spaces in between. Must be either this or None.
- substructure_search: Provide a SMARTS representation of a substructure to require in all candidates. If not using, must leave as None.
- number_substructure_search: Signifies how many occurrences of the substructure must appear in the candidate. If not a number, leave as None. If a substructure is given and this is left as None, the search will look for at least 1 or more occurrences.
- weights (optional): An array of weights signifying how to weigh each comparison value in the total rankings. This is [1,1,1,1,1] by default. These correspond to:
  - Structural Similarity
  - Molecular Weight Similarity
  - Thermophysical Similarity
  - Evaluated Toxicity
  - Synthetic Accessibility Scoring

Example Inputs:

    6517, 30, [True,True,False,False,False], False, [Si], CCO, 1

This example tells the CRS to search for recommendations for PubChem CID 6517. It would create a report with 30 candidates and use the thermophysical properties of Melting Point and Boiling Point in comparison. The elements allowed in candidates are the listed default plus Silicon. Finally, all returned candidates will have at least one occurrence of the SMARTS 'CCO'.

    quinolin-8-ol, 30, [False,False,False,False,True], True, None, None, None, [2,1,1,1,0]

This example tells the CRS to search for recommendations for the IUPAC Name quinolin-8-ol. It would create a report with 10 candidates and use the thermophysical property of the Henry's Law Constant. Candidates will be allowed to have any elements in it. The final sorting will disregard SA Scoring in the calculation and give increased weightage to structural similarity.
""",
            formatter_class=RawDescriptionHelpFormatter,
        )

        # Create a group for the mutually exclusive options -o and -t
        output_group = self.parser.add_mutually_exclusive_group(required=False)

        self.parser.add_argument(
            "-w",
            "--webapp",
            dest="webapp",
            help="Run the CRS interactively through the web-app in your browser",
            action="store_true",
        )

        output_group.add_argument(
            "-o", "--output", dest="output", help="Path to output report to be stored"
        )
        output_group.add_argument(
            "-t",
            "--time",
            dest="time",
            help="Generate output file with the current time as its name",
            action="store_true",
        )

        self.parser.add_argument(
            "-i",
            "--input",
            dest="input",
            help="Path to input txt file for run",
        )

        self.parser.add_argument(
            "-m",
            "--models",
            dest="containers",
            nargs="+",
            help="Names of models to be used",
            default=[],
        )

        args = self.parser.parse_args()

        if not args.webapp:
            if not (args.input and (args.output or args.time)):
                self.parser.error("If not using webapp, input and output name required")

    def parseArgs(self):
        """
        Parse the command-line arguments.

        Returns:
            argparse.Namespace: Parsed command-line arguments.
        """
        return self.parser.parse_args()

    def processData(self, args):
        """
        Process the input data based on the provided command-line arguments.

        Args:
            args (argparse.Namespace): Parsed command-line arguments.
        """
        if args.webapp:
            rf()
        else:
            with open(args.input, "r") as file:
                batch_input = file.read()

            # Generate output filename based on current time if -t is used
            if args.time:
                current_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
                output_filename = f"output/{current_time}"
            else:
                output_filename = "output/" + args.output

            # Process the containers if any
            containers = args.containers
            print(f"Using models: {containers}")  # Example usage, adjust as needed

            br(batch_input, output_filename, True, containers=containers)


if __name__ == "__main__":
    crs = CRS()
    args = crs.parseArgs()
    crs.processData(args)
