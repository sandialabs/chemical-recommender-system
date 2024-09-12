# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause


class AppState:
    """
    A class to maintain the state of the application.

    This class initializes various attributes that are used to store the state of the application,
    including captured output, chemical identifiers, query values, results, and other parameters.

    Attributes:
    - captured_output: Captured output to display progress.
    - cidarr: A list of chemical identifiers (CIDs) used in post-processing.
    - queryval: The query value provided by the user.
    - tarray: An array of boolean values indicating which thermophysical properties are to be calculated.
      The order is:
        - Melting Point (MP)
        - Boiling Point (BP)
        - Log P (logP)
        - Vapor Pressure (VP)
        - Henry's Law Constant (Hlaw)
    - results: The results of the comparison function.
    - query_models: Extra models used in the query.
    - finnum: The number of candidate results to be returned by the run.
    - params: Parameters for the comparison function.
    - subfailed: A flag indicating if substructure searching failed.
    - weights: An array of weights corresponding to the importance of each property in tarray.
    - batch_status: The status of the batch processing ("false", "loading", "true").
    - batch_process: The process running the batch processing.
    - search_status: The status of the search processing ("false", "loading", "true").
    - search_process: The process running the search processing.
    - result_queue: A queue to store the results from the search process.
    """

    def __init__(self, manager):
        self.captured_output = None
        self.cidarr = []
        self.queryval = ""
        self.tarray = []
        self.results = ""
        self.query_models = []
        self.finnum = None
        self.params = None
        self.subfailed = None
        self.weights = [1, 1, 1, 1, 1]
        self.batch_status = manager.Value("c", "false")
        self.batch_process = None
        self.search_status = manager.Value("c", "false")
        self.search_process = None
        self.result_queue = None

    def reset_vals(self):
        self.cidarr = []
        self.queryval = ""
        self.tarray = []
        self.results = ""
        self.query_models = []
        self.finnum = None
        self.params = None
        self.subfailed = None
        self.weights = [1, 1, 1, 1, 1]
        self.batch_process = None
        self.search_process = None
        self.search_status.value = "false"
        self.batch_status.value = "false"
        self.result_queue = None

        with open("logs/comparison-root.log", "w") as file:
            pass

        with open("logs/comparison.log", "w") as file:
            pass
