# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

from flask import Blueprint, render_template, request, current_app
import pandas, io, sys, shutil
import uuid
import os
from Comparison.utils.gen import parseQuery as pq
from Comparison.utils.graph import graphResults as gr
from Comparison.utils.pdf import combinePDFs as cp
from Comparison.utils.post import (
    createCheckstring as cc,
    processDict as pd,
    produceCSV as pc,
    formatResults as fr,
)
import traceback
import logging

SINGLE_JOBS_DIR = "src/Comparison/LocalIO/single_jobs"

results_bp = Blueprint("results", __name__)


@results_bp.route("/results", methods=["GET", "POST"])
def resultsFunc():
    logger = logging.getLogger(__name__)
    try:
        logger.info("Entered resultsFunc route.")
        # Defensive: try to get state, else use defaults
        try:
            state = current_app.state
        except Exception as e:
            logger.error(f"Could not get current_app.state: {e}")
            state = None

        # Defensive: get job_id
        job_id = request.args.get("job_id")
        if not job_id and state and hasattr(state, "single_job_id"):
            job_id = getattr(state, "single_job_id", None)
        logger.info(f"Using job_id: {job_id}")
        if not job_id:
            logger.error("No job_id provided or found in state.")
            return render_template(
                "results.html",
                error_message="No job ID provided. Cannot display results.",
                finalresults=[],
                tarray=[],
                weights=[1,1,1,1,1],
                results_metadata=[],
                job_id=None,
            ), 400

        job_dir = os.path.join(SINGLE_JOBS_DIR, job_id)
        logger.info(f"Job directory: {job_dir}")
        os.makedirs(job_dir, exist_ok=True)

        # Defensive: get weights
        weights = [1,1,1,1,1]
        if state and hasattr(state, "weights") and isinstance(state.weights, list) and len(state.weights) == 5:
            weights = state.weights
        elif state and hasattr(state, "weights"):
            try:
                weights = list(state.weights)
                while len(weights) < 5:
                    weights.append(1)
                weights = weights[:5]
            except Exception:
                weights = [1,1,1,1,1]

        # Handle POST request to update weights
        if request and request.method == "POST":
            logger.info("Processing POST request for weights.")
            for i in range(5):
                val = request.form.get(f"weight{i+1}")
                try:
                    weights[i] = float(val)
                except Exception:
                    weights[i] = 1

            if state:
                state.weights = weights

        # Defensive: get tarray
        tarray = []
        if state and hasattr(state, "tarray"):
            tarray = state.tarray

        # Defensive: get queryval
        queryval = ""
        if state and hasattr(state, "queryval"):
            queryval = state.queryval
            
        # Defensive: get results
        results = []
        results_metadata = []
        qcid = None
        opera_failed = False
        
        if state and hasattr(state, "results") and state.results:
            results = state.results
            logger.info(f"Found {len(results)} results in state.results")
        else:
            logger.warning("No results found in state, checking result queue")
            # Check if there are results in the queue that haven't been processed yet
            if (state and hasattr(state, "result_queue") and 
                state.result_queue is not None and not state.result_queue.empty()):
                logger.info("Found results in result_queue, processing...")
                result = state.result_queue.get()
                logger.info(f"Retrieved result from queue: status={result.get('status')}, job_id={result.get('job_id')}")
                
                if result["status"] == "true":
                    logger.info("Setting state values from successful result")
                    state.queryval = result["queryval"]
                    state.results = result["results"]
                    state.query_models = result["query_models"]
                    state.subfailed = result["subfailed"]
                    state.opera_failed = result.get("opera_failed", False)
                    state.search_status = "true"
                    
                    # Use the results we just retrieved
                    results = state.results
                    queryval = state.queryval
                    opera_failed = state.opera_failed
                    logger.info(f"Updated state - queryval: '{state.queryval}', results count: {len(state.results) if state.results else 0}, opera_failed: {state.opera_failed}")
                else:
                    logger.warning(f"Result status was not 'true'. Error: {result.get('error')}")
                    
        # Fallback: If no results in state or queue, try to read from CSV file directly
        if not results:
            logger.warning("No results found in state or queue, attempting to read from CSV file")
            try:
                import pandas as pd
                csv_path = "src/Comparison/LocalIO/data.csv"
                if os.path.exists(csv_path):
                    df = pd.read_csv(csv_path)
                    logger.info(f"Read CSV file with shape: {df.shape}")
                    if not df.empty:
                        # Convert DataFrame to list of lists (to match expected format)
                        results = df.values.tolist()
                        logger.info(f"Converted CSV to {len(results)} rows")
                        # Also update state if possible
                        if state:
                            state.results = results
                else:
                    logger.error(f"CSV file not found at {csv_path}")
            except Exception as csv_error:
                logger.error(f"Error reading CSV fallback: {csv_error}")
        
        logger.info(f"resultsFunc: Received results of type {type(results)} with length {len(results) if hasattr(results, '__len__') else 'N/A'}")
        logger.debug("=== RESULTS PROCESSING START ===")
        logger.debug(f"Job ID: {job_id}")
        logger.debug(f"Results type: {type(results)}")
        logger.debug(f"Results length: {len(results) if hasattr(results, '__len__') else 'N/A'}")
        
        if isinstance(results, list) and results and isinstance(results[0], (dict, list)):
            logger.debug(f"First result: {results[0]}")
            logger.info(f"resultsFunc: First result: {results[0]}")
            if isinstance(results[0], dict):
                logger.debug(f"First result keys: {list(results[0].keys())}")
                logger.info(f"resultsFunc: First result keys: {list(results[0].keys())}")
            elif isinstance(results[0], list):
                logger.debug(f"First result length: {len(results[0])}")
                logger.info(f"resultsFunc: First result length: {len(results[0])}")
        else:
            logger.debug("Results is empty or invalid format")
            logger.warning(f"Results is empty or invalid format: {results}")

        # Defensive: get finnum
        finnum = None
        if state and hasattr(state, "finnum"):
            finnum = state.finnum
        logger.info(f"resultsFunc: finnum = {finnum}")

        # Defensive: get tarray
        tarray = []
        if state and hasattr(state, "tarray"):
            tarray = state.tarray
        logger.info(f"resultsFunc: tarray = {tarray}")

        # Defensive: get weights
        weights = [1,1,1,1,1]
        if state and hasattr(state, "weights") and isinstance(state.weights, list) and len(state.weights) == 5:
            weights = state.weights
        elif state and hasattr(state, "weights"):
            try:
                weights = list(state.weights)
                while len(weights) < 5:
                    weights.append(1)
                weights = weights[:5]
            except Exception:
                weights = [1,1,1,1,1]
        logger.info(f"resultsFunc: weights = {weights}")

        # Defensive: get params
        params = None
        if state and hasattr(state, "params"):
            params = state.params

        # Defensive: get subfailed
        subfailed = None
        if state and hasattr(state, "subfailed"):
            subfailed = state.subfailed

        # Defensive: get opera_failed  
        opera_failed = False
        if state and hasattr(state, "opera_failed"):
            opera_failed = state.opera_failed

        logger.info("Creating checkstring from tarray.")
        checkstring = cc(tarray) if tarray else []

        logger.info("Parsing query value for CID.")
        logger.debug("Starting CID parsing process")
        try:
            qcid = (pq(queryval))[1]
            logger.debug(f"Parsed CID: {qcid}")
        except Exception:
            qcid = -1
        if qcid is None:
            qcid = -1

        logger.info("Processing results data.")
        finalresults = []
        fullcids = []
        try:
            # The results from SingleRun are already processed final_results, not raw comp_dict
            if isinstance(results, list) and results:
                logger.info(f"Processing {len(results)} results")
                finalresults = results
                
                # Extract CIDs from final results - include all results (query and candidates)
                fullcids = []
                for i, row in enumerate(results):
                    if row:  # Process all non-empty rows
                        # Extract CID from the result row
                        if isinstance(row, list) and len(row) > 0:
                            cid_raw = str(row[0])
                            
                            # Convert float string to int if needed, but preserve -1 for query
                            try:
                                if cid_raw == "-1":
                                    fullcids.append("-1")  # Keep query CID as is
                                elif '.' in cid_raw:
                                    cid_clean = str(int(float(cid_raw)))
                                    fullcids.append(cid_clean)
                                else:
                                    cid_clean = cid_raw.strip()
                                    fullcids.append(cid_clean)
                            except (ValueError, TypeError):
                                # If conversion fails, try to extract numeric part
                                import re
                                numeric_part = re.search(r'-?\d+', cid_raw)  # Include negative numbers
                                if numeric_part:
                                    fullcids.append(numeric_part.group())
                                else:
                                    fullcids.append(cid_raw.strip())
                        elif isinstance(row, dict) and 'cid' in row:
                            # Handle dictionary format
                            fullcids.append(str(row['cid']))
                            
                logger.info(f"Extracted {len(finalresults)} final results and {len(fullcids)} CIDs")
                logger.debug(f"Sample CIDs extracted: {fullcids[:3] if fullcids else 'None'}")
                
                # Debug: Print first few results to understand structure
                if finalresults:
                    logger.info(f"Sample results structure: {finalresults[:2]}")
            else:
                logger.warning(f"Results data is not in expected format: type={type(results)}, len={len(results) if hasattr(results, '__len__') else 'N/A'}")
        except Exception as e:
            logger.error(f"Error processing results data: {e}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            finalresults = []
            fullcids = []
        logger.debug("Results processing completed")
        logger.debug(f"Final results summary - finalresults: {len(finalresults) if finalresults else 0} items")
        # Create required directories
        for dir_path in ["src/Comparison/LocalIO", "src/App/static/LocalIO"]:
            os.makedirs(dir_path, exist_ok=True)

        logger.info("Copying Thermout.csv to static directory.")
        try:
            thermout_path = "src/Comparison/LocalIO/Thermout.csv"
            if os.path.exists(thermout_path):
                shutil.copy(
                    thermout_path,
                    "src/App/static/LocalIO/Property-Predictions.csv",
                )
            else:
                logger.warning(f"Thermout.csv not found at {thermout_path}")
                # Create empty file if missing
                with open("src/App/static/LocalIO/Property-Predictions.csv", "w") as f:
                    f.write("No property predictions available\n")
        except Exception as e:
            logger.error(f"Error copying Thermout.csv: {e}")

        # Generate graphs
        logger.info("Generating graphs.")
        try:
            # Verify data.csv exists before generating graphs
            if os.path.exists("src/Comparison/LocalIO/data.csv"):
                if qcid != -1:
                    gr(qcid, logger=logger)
                else:
                    gr(qcid, queryval, logger=logger)
                
                # Ensure graph is copied to the static directory
                if os.path.exists("src/Comparison/LocalIO/graph.png"):
                    shutil.copy(
                        "src/Comparison/LocalIO/graph.png",
                        "src/App/static/LocalIO/graph.png"
                    )
            else:
                logger.error("Cannot generate graph: data.csv not found")
        except Exception as e:
            logger.error(f"Error generating graphs: {e}")
            # Create a placeholder graph if we can't generate one
            try:
                import matplotlib.pyplot as plt
                plt.figure(figsize=(10, 6))
                plt.text(0.5, 0.5, "No graph data available", 
                        horizontalalignment='center', verticalalignment='center')
                plt.axis('off')
                plt.savefig("src/App/static/LocalIO/graph.png")
                logger.info("Created placeholder graph")
            except Exception as graph_e:
                logger.error(f"Could not create placeholder graph: {graph_e}")

        logger.info("Final results are ready for display.")
        # finalresults and fullcids are already available from the processing above

        logger.info("Formatting final results and updating CID array.")
        try:
            logger.debug("Formatting final results before processing")
            logger.debug(f"Pre-formatting results count: {len(finalresults) if finalresults else 0}")
            finalresults, cidarr = fr(finalresults, qcid)
            logger.debug("Final Results After formatting completed")
            logger.debug(f"Post-formatting results count: {len(finalresults) if finalresults else 0}")
            if state:
                state.cidarr = cidarr
                
            # Convert list-of-lists to list-of-objects for template compatibility
            if finalresults and isinstance(finalresults, list):
                formatted_results = []
                for row in finalresults:
                    if isinstance(row, list) and len(row) >= 7:
                        # Map list indices to named keys based on the data structure
                        # From debug output: [CID, ?, molecular_weight_sim, thermo_sim, ?, ?, structural_sim, ...]
                        import math
                        
                        def safe_float(val, default=0):
                            if val is None or val == "" or (isinstance(val, str) and val.lower() in ['nan', 'none']):
                                return default
                            try:
                                float_val = float(val)
                                return default if math.isnan(float_val) else float_val
                            except (ValueError, TypeError):
                                return default
                        
                        result_obj = {
                            "cid": row[0],                          # Column 0: CID (with formatting)
                            "final_score": safe_float(row[1] if len(row) > 1 else 0),        # Column 1: Final combined similarity score
                            "structural_similarity": safe_float(row[2] if len(row) > 2 else 0),  # Column 2: Structural similarity 
                            "mw_similarity": safe_float(row[3] if len(row) > 3 else 0),      # Column 3: MW similarity
                            "thermo_similarity": safe_float(row[4] if len(row) > 4 else 0),  # Column 4: Thermo similarity
                            "toxicity": safe_float(row[5] if len(row) > 5 else 0),          # Column 5: Toxicity
                            "sa_score": safe_float(row[6] if len(row) > 6 else 0),          # Column 6: SA Score
                        }
                        
                        # For CID formatting - ensure integers don't show as floats
                        if isinstance(row[0], str) and row[0] != "Query":
                            # Convert float CID to integer if needed
                            try:
                                if '.' in row[0] and row[0] != "Query":
                                    result_obj["cid"] = str(int(float(row[0])))
                            except (ValueError, TypeError):
                                pass  # Keep original if conversion fails
                        
                        # Add any additional columns from checkstring (property predictions)
                        for i, col in enumerate(checkstring):
                            col_index = 7 + i  # Properties start after the main similarity columns (CID + 6 values)
                            if len(row) > col_index:
                                val = row[col_index]
                                # For property predictions, keep original values but handle nan
                                if val is None or (isinstance(val, str) and val.lower() in ['nan', 'none']):
                                    result_obj[col] = ""
                                else:
                                    result_obj[col] = val
                            else:
                                result_obj[col] = ""
                                
                        formatted_results.append(result_obj)
                    else:
                        logger.warning(f"Skipping malformed result row: {row}")
                        
                finalresults = formatted_results
                
                # Sort by final similarity score in descending order, but keep the query (if any) at the top
                if finalresults:
                    query_results = []
                    non_query_results = []
                    
                    for result in finalresults:
                        # Check if this is the query compound (CID matches qcid or final_score is 0/None)
                        if (result.get("cid") == str(qcid) or 
                            result.get("final_score") == 0 or 
                            result.get("final_score") is None or
                            result.get("cid") == "Query"):
                            query_results.append(result)
                        else:
                            non_query_results.append(result)
                    
                    # Sort non-query results by final similarity score in descending order
                    non_query_results.sort(key=lambda x: x.get("final_score", 0), reverse=True)
                    
                    # Combine query results (at the top) with sorted non-query results
                    finalresults = query_results + non_query_results
                
                logger.info(f"Converted to {len(finalresults)} formatted result objects")
            else:
                finalresults = []
        except Exception as e:
            logger.error(f"Error in formatResults: {e}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            logger.error(f"finalresults at time of error: {finalresults}")
            logger.error(f"results at time of error: {results}")
            finalresults = []

        # Save CSV files after formatting is complete
        logger.info("Saving final results to CSV.")
        try:
            if finalresults:
                # For graph generation, use list format if available, otherwise convert
                if hasattr(results, '__iter__') and results and isinstance(results[0], list):
                    # Use original list format for graph
                    dfcsv_graph = pandas.DataFrame(results)
                    comp_data_path = "src/Comparison/LocalIO/data.csv"
                    dfcsv_graph.to_csv(comp_data_path, index=False)
                else:
                    logger.warning("Original list format not available for graph CSV")

                # For PDF generation, use formatted results with proper column names  
                if isinstance(finalresults[0], dict):
                    # Debug: Log the data being used for PDF
                    logger.info(f"DEBUG: First finalresult for PDF: {finalresults[0] if finalresults else 'None'}")
                    
                    # Create proper column mapping for PDF - match the UI table order exactly
                    pdf_data = []
                    for result in finalresults:
                        def round_if_numeric(val, decimals=4):
                            """Round numeric values to specified decimal places"""
                            if val is None or val == "" or (isinstance(val, str) and val.lower() in ['nan', 'none']):
                                return val
                            try:
                                if isinstance(val, (int, float)):
                                    return round(float(val), decimals)
                                elif isinstance(val, str) and val.replace('.', '').replace('-', '').isdigit():
                                    return round(float(val), decimals)
                                else:
                                    return val
                            except (ValueError, TypeError):
                                return val
                        
                        pdf_row = [
                            result.get('cid', ''),
                            round_if_numeric(result.get('final_score', '')),
                            round_if_numeric(result.get('structural_similarity', '')),  # Structural comes before MW to match UI table
                            round_if_numeric(result.get('mw_similarity', '')),
                            round_if_numeric(result.get('thermo_similarity', '')), 
                            round_if_numeric(result.get('toxicity', '')),
                            round_if_numeric(result.get('sa_score', '')),
                        ]
                        # Add any additional property columns with rounding
                        for col in checkstring:
                            pdf_row.append(round_if_numeric(result.get(col, '')))
                        pdf_data.append(pdf_row)
                    
                    logger.info(f"DEBUG: First PDF row: {pdf_data[0] if pdf_data else 'None'}")
                    
                    pdf_column_names = [
                        "PubChem CID",
                        "Final Similarity Score", 
                        "Structural Similarity",  # Match UI table order
                        "Molecular Weight Similarity",
                        "Thermophysical Similarity",
                        "Predicted Toxicity",
                        "Synthetic Availability Score",
                    ]
                    # Add property column names
                    for col in checkstring:
                        label_map = {
                            "MP_pred": "Predicted Melting Point (°C)",
                            "BP_pred": "Predicted Boiling Point (°C)",
                            "LogP_pred": "Predicted LogP",
                            "LogHL_pred": "Predicted Henry's Law Constant",
                            "LogVP_pred": "Predicted Vapor Pressure",
                            "LogBCF_pred": "Predicted BCF",
                            "CATMoS_EPA_pred": "CATMoS EPA",
                            "CATMoS_LD50_pred": "LD50",
                        }
                        pdf_column_names.append(label_map.get(col, col))
                    
                    dfcsv_pdf = pandas.DataFrame(pdf_data, columns=pdf_column_names)
                    
                    # Save CSV for PDF with proper headers
                    csv_path = os.path.join(job_dir, "results.csv")
                    dfcsv_pdf.to_csv(csv_path, index=False)
                    
                    # Copy to static directory for PDF generation
                    static_results_csv = os.path.join("src", "App", "static", "LocalIO", "results.csv")
                    os.makedirs(os.path.dirname(static_results_csv), exist_ok=True)
                    shutil.copy(csv_path, static_results_csv)
                    
                    logger.info(f"CSV saved to {csv_path} and {static_results_csv} (with proper headers for PDF)")
                else:
                    logger.warning("Final results not in expected object format")
            else:
                logger.warning("No final results to save to CSV")
        except Exception as e:
            logger.error(f"Error saving results.csv: {e}")

        logger.info("Combining PDFs for final report.")
        try:
            cp(
                fullcids,
                queryval,
                qcid,
                params,
                weights,
                True,
                subfailed,
            )
        except Exception as e:
            logger.error(f"Error in combinePDFs: {e}")

        logger.info("Building results metadata.")
        # Debug: Log the thermophysical values being sent to the template
        if finalresults:
            thermo_values = []
            for i, result in enumerate(finalresults):
                if isinstance(result, dict):
                    thermo_val = result.get('thermo_similarity', 'N/A')
                    thermo_values.append(f"Row {i}: {thermo_val}")
                elif isinstance(result, list) and len(result) > 3:
                    thermo_val = result[3]
                    thermo_values.append(f"Row {i}: {thermo_val}")
            logger.info(f"DEBUG: Thermophysical values being sent to table: {thermo_values}")
        
        results_metadata = [
            {"key": "cid", "label": "PubChem CID", "type": "string"},
            {"key": "final_score", "label": "Final Similarity Score", "type": "number"},
            {"key": "structural_similarity", "label": "Structural Similarity", "type": "number"},
            {"key": "mw_similarity", "label": "Molecular Weight Similarity", "type": "number"},
            {"key": "thermo_similarity", "label": "Thermophysical Similarity", "type": "number"},
            {"key": "toxicity", "label": "Predicted Toxicity", "type": "number"},
            {"key": "sa_score", "label": "Synthetic Availability Score", "type": "number"},
        ]
        for col in checkstring:
            label_map = {
                "MP_pred": "Predicted Melting Point (°C)",
                "BP_pred": "Predicted Boiling Point (°C)",
                "LogP_pred": "Predicted LogP",
                "LogHL_pred": "Predicted Henry's Law Constant",
                "LogVP_pred": "Predicted Vapor Pressure",
                "LogBCF_pred": "Predicted BCF",
                "CATMoS_EPA_pred": "CATMoS EPA",
                "CATMoS_LD50_pred": "LD50",
            }
            results_metadata.append({
                "key": col,
                "label": label_map.get(col, col),
                "type": "string"
            })

        logger.info("Linking/copying job_dir to static/single_jobs/<job_id>.")
        static_jobs_dir = os.path.join(current_app.root_path, "static", "single_jobs")
        os.makedirs(static_jobs_dir, exist_ok=True)
        static_job_dir = os.path.join(static_jobs_dir, job_id)
        if not os.path.exists(static_job_dir):
            try:
                os.symlink(os.path.abspath(job_dir), static_job_dir)
            except Exception:
                try:
                    shutil.copytree(job_dir, static_job_dir)
                except Exception as e:
                    logger.error(f"Error linking/copying job_dir: {e}")

        logger.info("Rendering results.html template.")
        logger.info(f"resultsFunc: finalresults type: {type(finalresults)}, length: {len(finalresults) if hasattr(finalresults, '__len__') else 'N/A'}")
        
        # Final debugging before template render
        logger.debug("=== FINAL TEMPLATE DATA ===")
        logger.debug(f"finalresults type: {type(finalresults)}")
        logger.debug(f"finalresults length: {len(finalresults) if hasattr(finalresults, '__len__') else 'No length'}")
        logger.debug(f"results_metadata type: {type(results_metadata)}")
        logger.debug(f"results_metadata length: {len(results_metadata) if hasattr(results_metadata, '__len__') else 'No length'}")
        logger.debug(f"queryval: {queryval}")
        logger.debug(f"qcid: {qcid}")
        logger.debug("=== RENDERING TEMPLATE ===")
        
        return render_template(
            "results.html",
            finalresults=finalresults,
            tarray=tarray,
            weights=weights,
            results_metadata=results_metadata,
            job_id=job_id,
            qcid=qcid,
            state=state,
            opera_failed=opera_failed,  # Pass OPERA failure status to template
        )
    except Exception as e:
        tb = traceback.format_exc()
        logger.error(f"Exception in resultsFunc: {e}\n{tb}")
        return render_template(
            "results.html",
            error_message=f"An error occurred while generating results: {e}\n\n{tb}",
            finalresults=[],
            tarray=[],
            weights=[1,1,1,1,1],
            results_metadata=[],
            job_id=None,
            opera_failed=False,  # Default to false on error
        ), 500
