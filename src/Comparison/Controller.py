# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pandas, os
from pypdf import PdfWriter
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

# Suppress watchdog debug logging immediately
import logging
logging.getLogger("watchdog").setLevel(logging.WARNING)
logging.getLogger("watchdog.observers").setLevel(logging.WARNING)
logging.getLogger("watchdog.observers.inotify_buffer").setLevel(logging.WARNING)

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

from utils.progress_logger import progress_context
import traceback


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


def monitor_log(log_file):
    event_handler = LogHandler(log_file)
    observer = Observer()
    observer.schedule(event_handler, path=os.path.dirname(log_file), recursive=False)
    observer.start()
    return observer


def SingleRun(
    queryinput, finnum, tarray, incEle, smarts, smarts_num, weights, containers, job_id="default", progress_queue=None, disallow_isotopes=False
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
    - job_id: Unique identifier for this job (for logging and progress tracking).
    - progress_queue: Optional queue for real-time progress updates to UI.

    Returns:
    - Tuple containing (results, query_models, subfailed) or None if there was an error.
    """
    log_file_path = "logs/comparison-root.log"
    observer = monitor_log(log_file_path)

    with progress_context(job_id, progress_queue) as logger:

        logger.step("Starting comparison function", "Running")
        try:
            logger.step("Validating inputs and parameters", "Running")
            
            # === TECHNICAL VALIDATION SECTION ===
            # Validate tarray structure and indices
            if not isinstance(tarray, (list, tuple)) or len(tarray) != 5:
                raise ValueError(f"tarray must be a list/tuple of exactly 5 boolean values, got: {tarray}")
            if not all(isinstance(x, bool) for x in tarray):
                raise ValueError(f"All tarray elements must be boolean, got: {[type(x).__name__ for x in tarray]}")
            
            # Validate weights array structure  
            expected_weights_length = 5 + (len(containers) if containers else 0)
            if not isinstance(weights, (list, tuple)) or len(weights) != expected_weights_length:
                raise ValueError(f"weights array must have {expected_weights_length} elements (5 base + {len(containers) if containers else 0} containers), got {len(weights) if weights else 'None'}: {weights}")
            if not all(isinstance(x, (int, float)) and x >= 0 for x in weights):
                raise ValueError(f"All weight values must be non-negative numbers, got: {weights}")
                
            # Validate finnum parameter
            if not isinstance(finnum, int) or finnum <= 0:
                raise ValueError(f"finnum must be a positive integer, got: {finnum}")
            if finnum > 1000:  # Reasonable upper limit
                logger.warning(f"finnum={finnum} is very large, this may cause performance issues")
                
            # Validate smarts_num consistency
            # Convert smarts_num to int if possible
            if smarts and smarts_num is not None:
                try:
                    smarts_num = int(smarts_num)
                except Exception:
                    raise ValueError(f"smarts_num must be a non-negative integer when smarts is provided, got: {smarts_num}")
                if smarts_num < 0:
                    raise ValueError(f"smarts_num must be a non-negative integer when smarts is provided, got: {smarts_num}")
            if smarts_num and not smarts:
                logger.warning("smarts_num provided without smarts pattern, will be ignored")
                
            # Validate container names
            if containers and not isinstance(containers, (list, tuple)):
                raise ValueError(f"containers must be a list or tuple, got: {type(containers)}")
            if containers and not all(isinstance(x, str) for x in containers):
                raise ValueError(f"All container names must be strings, got: {[type(x).__name__ for x in containers]}")
                
            logger.info(f"Input validation passed - query={queryinput}, finnum={finnum}, tarray={tarray}, weights={weights}, containers={containers}")
            # === END VALIDATION SECTION ===
            
            logger.step("Parsing query and preparing search parameters", "Running")
            logger.info(f"SingleRun called with: query={queryinput}, finnum={finnum}, job_id={job_id}")
            
            # Begin single run by calling the overarching Comparison Function
            logger.progress("Running", "Finding similar compounds...")
            logger.info("About to call cf() - comparison function")
            
            data = cf(
                queryinput,
                finnum,
                tarray,
                incEle,
                smarts,
                smarts_num,
                containers=containers,
                job_id=job_id,
                disallow_isotopes=disallow_isotopes,
            )
            
            # === TECHNICAL VALIDATION OF CF RETURN ===
            if not isinstance(data, (list, tuple)):
                raise ValueError(f"comparisonFunction must return [comp_dict, query_models, subfailed, opera_failed], got: {type(data)}")
            
            # Handle backward compatibility - if only 3 values returned, assume opera_failed=False
            if len(data) == 3:
                logger.warning("comparisonFunction returned 3 values instead of 4 - assuming opera_failed=False for backward compatibility")
                data = list(data) + [False]  # Add default opera_failed=False
            elif len(data) != 4:
                raise ValueError(f"comparisonFunction must return [comp_dict, query_models, subfailed, opera_failed], got: {len(data)} values")
            
            logger.info(f"cf() completed successfully, data type: {type(data)}, data length: {len(data) if data else 'None'}")
            logger.progress("Processing", "Analyzing and ranking results...")

            # Initialize data values from comparison results
            results = data[0]  # comp_dict
            query_models = data[1]  # query_models array
            subfailed = data[2]  # subfailed boolean
            opera_failed = data[3]  # opera_failed boolean
            
            # Validate the structure of returned data
            if not isinstance(results, dict):
                raise ValueError(f"results (comp_dict) must be a dictionary, got: {type(results)}")
            if query_models is not None and not isinstance(query_models, (list, tuple)):
                raise ValueError(f"query_models must be a list/tuple or None, got: {type(query_models)}")
            if not isinstance(subfailed, bool):
                raise ValueError(f"subfailed must be a boolean, got: {type(subfailed)}")
                
            # Validate comp_dict structure - each value should be [score, structural, molecular, thermal, toxic, SA, ...containers]
            for cid, values in results.items():
                if not isinstance(values, (list, tuple)):
                    raise ValueError(f"comp_dict values must be lists/tuples, CID {cid} has: {type(values)}")
                expected_length = 6 + (len(containers) if containers else 0)
                if len(values) != expected_length:
                    raise ValueError(f"comp_dict values must have {expected_length} elements, CID {cid} has {len(values)}: {values}")

            logger.info(f"cf() returned results: {len(results) if results else 0} items, type: {type(results)}")
            logger.info(f"query_models: {len(query_models) if query_models else 0} items, type: {type(query_models)}")
            logger.info(f"subfailed: {subfailed}")
            logger.info(f"opera_failed: {opera_failed}")
            
            # Log OPERA failure warning if needed
            if opera_failed:
                logger.warning("OPERA property predictions failed - thermal/toxicity comparisons use neutral values")
                logger.progress("Warning", "Property calculations failed - analysis continues with limited accuracy")
            # === END CF VALIDATION ===
            
            if not results or len(results) == 0:
                logger.error("No candidates found. Please check your query or try different parameters.")
                logger.progress("Error", "No candidates found. Please check your query or try different parameters.")
                return None

            # Check if only the query itself is present (comp_dict is keyed by CID strings)
            candidate_cids = list(results.keys())
            logger.info(f"Found {len(candidate_cids)} candidate CIDs: {candidate_cids[:5]}..." if len(candidate_cids) > 5 else f"Found {len(candidate_cids)} candidate CIDs: {candidate_cids}")
            
            if len(candidate_cids) == 1 and str(candidate_cids[0]) == str(queryinput):
                logger.error("Only the query itself was found. No candidates generated.")
                logger.progress("Error", "Only the query itself was found. No candidates generated.")
                return None

            logger.progress("Processing", "Finding PubChem CID for the query.")
            qcid = (pq(queryinput))[1]
            if qcid is None:
                qcid = -1
            logger.info(f"Query CID resolved to: {qcid}")

            logger.progress("Processing", "Determining which properties were calculated.")
            checkstring = cc(tarray)
            logger.info(f"Property check string: {checkstring}")

            logger.progress("Processing", "Calculating similarity scores...")
            logger.info(f"About to call processDict with {len(results)} results, finnum={finnum}")
            processed_dict = pd(results, weights, qcid, finnum, query_models, containers)
            logger.info(f"processDict returned {len(processed_dict) if processed_dict else 0} candidates.")
            logger.info(f"processed_dict type: {type(processed_dict)}, sample keys: {list(processed_dict.keys())[:3] if processed_dict else 'None'}")

            if not processed_dict or len(processed_dict) == 0:
                logger.error("No processed candidates available after ranking. Please check your input.")
                logger.progress("Error", "No processed candidates available after ranking.")
                return None

            logger.progress("Processing", "Saving results to CSV.")
            try:
                # First generate the final results from produceCSV
                logger.info(f"About to call produceCSV with processed_dict ({len(processed_dict)} items), checkstring={checkstring}")
                final_results, fullcids = pc(processed_dict, checkstring, tarray, containers, query_models)
                logger.info(f"produceCSV returned {len(fullcids) if fullcids else 0} CIDs and {len(final_results) if final_results else 0} formatted results.")
                logger.info(f"final_results type: {type(final_results)}, fullcids type: {type(fullcids)}")
                if final_results and len(final_results) > 0:
                    logger.info(f"Sample final_result: {final_results[0]}")
                
                # Create DataFrame from the final results list
                if final_results:
                    try:
                        # Define proper column names for the CSV based on the expected structure
                        # final_results is a list of lists from produceCSV
                        expected_columns = [
                            'CID',           # 0 - Compound ID
                            'Overall',       # 1 - val[0] - Overall weighted similarity score
                            'Fingerprint',   # 2 - val[1] - Fingerprint/structural (Tanimoto) similarity  
                            'Molecular',     # 3 - val[2] - Molecular weight similarity
                            'Thermophysical', # 4 - val[3] - Thermophysical similarity
                            'Toxicity',      # 5 - val[4] - Toxicity score
                            'Synthetic'      # 6 - val[5] - Synthetic availability
                        ]
                        
                        # Extend columns list if there are more columns in the data
                        if final_results:
                            num_cols = len(final_results[0])
                            while len(expected_columns) < num_cols:
                                expected_columns.append(f'Extra_{len(expected_columns)}')
                        
                        dfcsv = pandas.DataFrame(final_results, columns=expected_columns[:len(final_results[0]) if final_results else 8])
                        dfcsv.to_csv(r"src/Comparison/LocalIO/data.csv", index=False)
                        logger.info(f"Saved data.csv with shape {dfcsv.shape}, columns: {list(dfcsv.columns)}")
                    except Exception as csv_error:
                        logger.error(f"Error creating DataFrame from final_results: {csv_error}")
                        # Create a minimal DataFrame as fallback
                        fallback_data = [{"CID": "N/A", "Structural": 0, "Fingerprint": 0, "Thermophysical": 0, "Toxicity": 0, "Synthetic": 1}]
                        dfcsv = pandas.DataFrame(fallback_data)
                        dfcsv.to_csv(r"src/Comparison/LocalIO/data.csv", index=False)
                        logger.warning("Created fallback data.csv with proper structure")
                else:
                    logger.warning("No final results to save to CSV.")
                    # Create empty CSV with proper structure
                    fallback_data = [{"CID": "N/A", "Structural": 0, "Fingerprint": 0, "Thermophysical": 0, "Toxicity": 0, "Synthetic": 1}]
                    dfcsv = pandas.DataFrame(fallback_data)
                    dfcsv.to_csv(r"src/Comparison/LocalIO/data.csv", index=False)
            except Exception as e:
                logger.error(f"Error saving data.csv: {e}")
                # Ensure we have some CSV file for downstream processes
                try:
                    # Create a fallback CSV with the correct structure for graph generation
                    fallback_data = [{"CID": "N/A", "Structural": 0, "Fingerprint": 0, "Thermophysical": 0, "Toxicity": 0, "Synthetic": 1}]
                    dfcsv = pandas.DataFrame(fallback_data)
                    dfcsv.to_csv(r"src/Comparison/LocalIO/data.csv", index=False)
                    logger.info("Created emergency fallback data.csv with proper structure")
                except:
                    pass
                final_results, fullcids = [], []

            logger.progress("Processing", "Generating visualizations...")
            try:
                if qcid == -1:
                    gr(qcid, queryinput, logger=logger)
                else:
                    gr(qcid, logger=logger)
                logger.info("Graph generated successfully.")
            except Exception as e:
                logger.error(f"Error generating graph: {e}")

            logger.progress("Processing", "Creating PDF report...")
            try:
                params = [queryinput, finnum, tarray, incEle, smarts, smarts_num, disallow_isotopes]
                cp(
                    fullcids,
                    queryinput,
                    qcid,
                    params,
                    weights,
                    False,
                    subfailed,
                    containers=containers,
                    opera_failed=opera_failed,  # Pass OPERA failure status to PDF generation
                )
                logger.info("PDF report generated successfully.")
            except Exception as e:
                logger.error(f"Error generating PDF report: {e}")

            logger.progress("Completed", "Analysis complete!")
            logger.info(f"[DEBUG] About to return from SingleRun: final_results length={len(final_results) if final_results else 0}")
            logger.info(f"[DEBUG] final_results type: {type(final_results)}")
            if final_results and len(final_results) > 0:
                logger.info(f"[DEBUG] First result sample: {final_results[0] if len(final_results) > 0 else 'None'}")
            
            # Return the final formatted results for use by the caller (not the raw comp_dict)
            return (final_results, query_models, subfailed, opera_failed)
        except Exception as e:
            tb = traceback.format_exc()
            logger.error(f"Exception in SingleRun: {e}\n{tb}")
            # Optionally, re-raise if you want Flask to handle it as 500
            return None  # Return None on error
        finally:
            observer.stop()
            observer.join()
            # Don't send a final progress message here as it might interfere with batch processing


def BatchRun(batch_text, output=None, from_command=False, containers=[], job_id="default", progress_queue=None):
    """
    Perform a batch run of the comparison function for multiple queries.

    Inputs:
    - batch_text: The input text containing multiple queries.
    - output: The output file path for the final PDF report (default is None).
    - from_command: A flag indicating whether the input is from a command line text file (default is False).
    - containers: A list of container names for additional models (default is None).
    - job_id: Unique identifier for this job (for logging and progress tracking).
    - progress_queue: Optional queue for real-time progress updates to UI.

    Returns:
    - None (generates and saves the necessary outputs).
    """
    with progress_context(job_id, progress_queue) as logger:
        # Line by line parse the input in batch text, pull parameter values for a singular search
        # The starting page of all reports should be the CRS summary
        merger = PdfWriter()
        merger.append("src/App/static/assets/CRS1pagesum.pdf")

        logger.step("Parsing batch input", "Processing")

        # Parse the batch input text
        queries = pb(batch_text, from_command, containers, job_id)

        ran = False
        dfs = []

        # Perform a single run for each parsed query
        run_num = 0
        total_queries = len(queries)
        logger.info(f"Starting batch processing of {total_queries} queries")
        logger.progress("Initializing", f"Batch job started: {total_queries} queries to process")
        
        for query_params in queries:
            run_num = run_num + 1
            
            # More detailed progress message showing which query is being processed
            query_display = query_params['query']
            if len(query_display) > 50:
                query_display = query_display[:47] + "..."
            
            logger.progress("Processing", f"Query {run_num} of {total_queries}: {query_display}")
            
            try:
                # Use the main job_id (not sub_job_id) so progress updates show in the main stream
                # But still pass a unique identifier for logging
                sub_job_id = f"{job_id}_query_{run_num}"
                
                logger.progress("Processing", f"Starting analysis for query {run_num}: {query_display}")
                
                SingleRun_result = SingleRun(
                    query_params["query"],
                    query_params["finnum"],
                    query_params["tarray"],
                    query_params["incEle"],
                    query_params["smarts"],
                    query_params["smarts_num"],
                    query_params["weights"],
                    containers,
                    job_id,  # Use main job_id instead of sub_job_id for progress updates
                    progress_queue,
                    query_params.get("disallow_isotopes", False),
                )
                
                # Check for OPERA failure in this run
                if SingleRun_result and len(SingleRun_result) >= 4:
                    opera_failed = SingleRun_result[3]
                    if opera_failed:
                        logger.progress("Warning", f"Query {run_num}: Property predictions failed - results may have reduced accuracy")
                
                ran = True
                merger.append("src/App/static/LocalIO/report.pdf")
                dfs.append(pandas.DataFrame({"Run": [f"run {run_num}"]}))
                
                # Check if Thermout.csv indicates OPERA failure
                thermout_df = pandas.read_csv("src/Comparison/LocalIO/Thermout.csv")
                if 'OPERA_FAILED' in thermout_df.columns and thermout_df['OPERA_FAILED'].any():
                    logger.warning(f"Query {run_num}: OPERA failed - thermal/toxicity data is synthetic")
                    # Add a note to the dataframe
                    thermout_df['OPERA_STATUS'] = 'FAILED'
                else:
                    thermout_df['OPERA_STATUS'] = 'SUCCESS'
                    
                dfs.append(thermout_df)
                
                # Progress update after successful completion of this query
                logger.progress("Processing", f"Completed query {run_num} of {total_queries}: {query_display}")
                
            except Exception as e:
                logger.error(f"Error processing query {run_num} ({query_params['query']}): {str(e)}")
                logger.progress("Processing", f"Error in query {run_num} of {total_queries}: {query_display} - continuing with next query")
                # Continue with other queries even if one fails
                continue

        logger.progress("Processing", "Combining results and generating final report")

        # If output is provided, use it as the base path for all outputs (CSV, PDF, metadata)
        if output:
            base_path = os.path.splitext(output)[0]
            csv_path = base_path + ".csv"
            pdf_path = base_path + ".pdf"
            metadata_path = base_path + "-metadata.json"
        else:
            csv_path = "src/App/static/LocalIO/Combined-OPERA-Results.csv"
            pdf_path = "src/App/static/LocalIO/Batch-Report.pdf"
            metadata_path = "src/App/static/LocalIO/Combined-OPERA-Results-metadata.json"

        # Concatenate all DataFrames in the list into a single DataFrame
        if not dfs:
            logger.error("No data was collected from any queries - cannot generate combined results")
            logger.progress("Error", "Batch processing failed - no data was collected from any queries")
            # Create a minimal error DataFrame
            combined_df = pandas.DataFrame({
                "Error": ["No queries produced valid results"],
                "BatchStatus": ["FAILED"]
            })
        else:
            combined_df = pandas.concat(dfs, ignore_index=True)
            logger.info(f"Successfully combined data from {len(dfs)} dataframes")
            
        combined_df.to_csv(csv_path, index=False)

        logger.info(f"Combined results saved to {csv_path}")

        # Generate and save metadata for the batch results
        from Comparison.utils.gen import generate_metadata_from_dataframe
        metadata = generate_metadata_from_dataframe(combined_df)
        import json
        with open(metadata_path, "w") as f:
            json.dump(metadata, f)

        # Write resultant PDFs together and finally save to return path
        if ran:
            merger.write(pdf_path)
            logger.info(f"Batch report saved to {pdf_path}")
            logger.progress("Finished", f"All {total_queries} queries completed successfully! Results are ready for download.")
        else:
            logger.error("No queries were successfully processed")
            logger.progress("Error", "Batch processing failed - no queries were successfully processed")
            logger.warning("No queries were processed successfully")
        
        if output:
            # Already written above
            pass
        merger.close()


# Testing Runs:
# SingleRun(1923, 10, [False,True,False,True,False], False, True, None, None)
# SingleRun('C2(CC1(=CC=C(C=C1)N))(=CC=C(C=C2)N)', 10, [False,True,False,True,False], False, True, '[NH2]', 2)
# SingleRun("CCCCCCCCC2C(CCCCCCCC(=O)OCC1CO1)C=CC(CCCCCC)C2CCCCCCCC(=O)OCC3CO3", 10,  [True,True,True,True,True], False, True, None, None) # <- run that doesn't exist in PubChem
