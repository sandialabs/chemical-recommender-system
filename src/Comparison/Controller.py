# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pandas, time, os
from pypdf import PdfWriter
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

from Comparison.Comparison import comparisonFunction as cf

from Comparison.utils.gen import parseQuery as pq, parseBatchInput as pb
from Comparison.utils.graph import graphResults as gr
from Comparison.utils.norm import normalizeDict as nd
from Comparison.utils.csv import createCSVOut as cv
from Comparison.utils.pdf import combinePDFs as cp
from Comparison.utils.post import (
    createCheckstring as cc,
    processDict as pd,
    produceCSV as pc,
)


class LogHandler(FileSystemEventHandler):
    def __init__(self, log_file):
        self.log_file = log_file
        self.last_position = 0

    def on_modified(self, event):
        if event.src_path == self.log_file:
            # Check the current size of the file
            current_size = os.path.getsize(self.log_file)

            # If the file size is smaller than the last known position, reset the position
            if current_size < self.last_position:
                self.last_position = 0

            with open(self.log_file, "r") as file:
                # Move to the last read position
                file.seek(self.last_position)
                # Read new lines
                new_lines = file.readlines()
                # Update the last read position
                self.last_position = file.tell()
                # Print new lines to the terminal
                for line in new_lines:
                    print(line, end="")


def monitor_log(log_file):
    event_handler = LogHandler(log_file)
    observer = Observer()
    observer.schedule(event_handler, path=os.path.dirname(log_file), recursive=False)
    observer.start()
    return observer


def SingleRun(
    queryinput, finnum, tarray, incEle, smarts, smarts_num, weights, containers
):
    """
    Perform a single run of the comparison function and generate the necessary outputs.

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
    - incEle: A flag indicating whether to include all elements.
    - smarts: The SMARTS pattern for substructure searching.
    - smarts_num: The number of substructure matches required.
    - weights: A list of weights for the query.
    - containers: A list of container names for additional models.

    Returns:
    - None (generates and saves the necessary outputs).
    """
    log_file_path = "logs/comparison-root.log"
    observer = monitor_log(log_file_path)

    try:
        # Begin single run by calling the overarching Comparison Function
        data = cf(
            queryinput,
            finnum,
            tarray,
            incEle,
            smarts,
            smarts_num,
            containers=containers,
        )
    finally:
        observer.stop()
        observer.join()

    # Initialize data values from comparison results
    results = data[0]
    query_models = data[1]
    subfailed = data[2]

    # Find a CID for the query. If it doesn't exist in PubChem, label as -1
    qcid = (pq(queryinput))[1]
    if qcid is None:
        qcid = -1

    # Create an array of which OPERA properties were calculated from search parameters
    checkstring = cc(tarray)

    # Post-processing on dict, return normalized top n values as asked for from search parameters
    processed_dict = pd(results, weights, qcid, finnum, query_models, containers)

    # Store data information and use it to create graphs of data
    dfcsv = pandas.DataFrame(processed_dict)
    dfcsv.to_csv(r"src/Comparison/LocalIO/data.csv", index=False)

    if qcid == -1:
        gr(qcid, queryinput)
    else:
        gr(qcid)

    # Create and save a CSV of the results for user download, simultaneously get a list of resultant CIDs
    _, fullcids = pc(processed_dict, checkstring, tarray, containers, query_models)

    # Create the PDF report of the entire run, store it in LocalIO
    params = [queryinput, finnum, tarray, incEle, smarts, smarts_num]
    cp(
        fullcids,
        queryinput,
        qcid,
        params,
        weights,
        False,
        subfailed,
        containers=containers,
    )


def BatchRun(batch_text, output=None, from_command=False, containers=[]):
    """
    Perform a batch run of the comparison function for multiple queries.

    Inputs:
    - batch_text: The input text containing multiple queries.
    - output: The output file path for the final PDF report (default is None).
    - from_command: A flag indicating whether the input is from a command line text file (default is False).
    - containers: A list of container names for additional models (default is None).

    Returns:
    - None (generates and saves the necessary outputs).
    """
    # Line by line parse the input in batch text, pull parameter values for a singular search
    # The starting page of all reports should be the CRS summary
    merger = PdfWriter()
    merger.append("src/App/static/assets/CRS1pagesum.pdf")

    # Parse the batch input text
    queries = pb(batch_text, from_command, containers)

    ran = False
    dfs = []

    # Perform a single run for each parsed query
    run_num = 0
    for query_params in queries:
        run_num = run_num + 1
        SingleRun(
            query_params["query"],
            query_params["finnum"],
            query_params["tarray"],
            query_params["incEle"],
            query_params["smarts"],
            query_params["smarts_num"],
            query_params["weights"],
            containers,
        )
        ran = True
        merger.append("src/App/static/LocalIO/report.pdf")
        dfs.append(pandas.DataFrame({"Run": [f"run {run_num}"]}))
        dfs.append(pandas.read_csv("src/Comparison/LocalIO/Thermout.csv"))

    # Concatenate all DataFrames in the list into a single DataFrame
    combined_df = pandas.concat(dfs, ignore_index=True)
    # Save the combined DataFrame to a new CSV file
    combined_df.to_csv("src/App/static/LocalIO/Combined-OPERA-Results.csv", index=False)

    # Write resultant PDFs together and finally save to return path
    if ran:
        merger.write("src/App/static/LocalIO/Batch-Report.pdf")
    if output:
        merger.write(output + ".pdf")
        combined_df.to_csv(output + ".csv", index=False)
    merger.close()


# Testing Runs:
# SingleRun(1923, 10, [False,True,False,True,False], False, True, None, None)
# SingleRun('C2(CC1(=CC=C(C=C1)N))(=CC=C(C=C2)N)', 10, [False,True,False,True,False], False, True, '[NH2]', 2)
# SingleRun("CCCCCCCCC2C(CCCCCCCC(=O)OCC1CO1)C=CC(CCCCCC)C2CCCCCCCC(=O)OCC3CO3", 10,  [True,True,True,True,True], False, True, None, None) # <- run that doesn't exist in PubChem
