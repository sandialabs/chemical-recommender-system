# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import time
import numpy as np
import os
from pymilvus import MilvusClient
from utils.progress_logger import get_progress_logger
from flask import current_app

logger = None  # Will be replaced with job-specific logger

# Create necessary directories to prevent file save errors
def ensure_directories_exist(job_id="default"):
    """
    Create all required directories for the application to run properly.
    This prevents "Cannot save file into a non-existent directory" errors.
    """
    logger = get_progress_logger(job_id)
    directories = [
        "src/Comparison/LocalIO",
        "src/App/static/LocalIO",
        "src/App/static/jobs"
    ]
    
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        logger.info(f"Ensured directory exists: {directory}")


def connect_to_milvus_with_retry(job_id="default", max_retries=3, retry_delay=30):
    """
    Connect to Milvus with retry logic in case the container is restarting.
    
    Args:
        job_id: Unique identifier for this job (for logging purposes)
        max_retries: Maximum number of connection attempts
        retry_delay: Delay between retry attempts in seconds
        
    Returns:
        MilvusClient: Connected Milvus client
        
    Raises:
        Exception: If connection fails after all retries
    """
    logger = get_progress_logger(job_id)
    
    for attempt in range(max_retries):
        try:
            logger.info(f"Attempting to connect to Milvus (attempt {attempt + 1}/{max_retries})...")
            
            # Try localhost first
            try:
                client = MilvusClient(
                    uri="tcp://localhost:19530",
                    timeout=30,  # 30 second timeout
                    # Add connection pooling parameters if available
                )
                logger.info("Successfully connected to Milvus Client on tcp://localhost:19530")
                return client
            except Exception as e:
                logger.warning(f"Failed to connect to localhost:19530: {str(e)}")
                
                # Try standalone container
                try:
                    client = MilvusClient(
                        uri="tcp://standalone:19530",
                        timeout=30,  # 30 second timeout
                    )
                    logger.info("Successfully connected to Milvus Client on tcp://standalone:19530")
                    return client
                except Exception as e2:
                    logger.warning(f"Failed to connect to standalone:19530: {str(e2)}")
                    raise e2
                    
        except Exception as e:
            logger.error(f"Connection attempt {attempt + 1} failed: {str(e)}")
            
            if attempt < max_retries - 1:
                logger.info(f"Waiting {retry_delay} seconds before retry...")
                time.sleep(retry_delay)
                
                # Force garbage collection between connection attempts
                import gc
                gc.collect()
                logger.debug("Forced garbage collection between connection attempts")
            else:
                logger.error(f"All {max_retries} connection attempts failed")
                raise Exception(f"Failed to connect to Milvus after {max_retries} attempts: {str(e)}")


def search_partitions_with_retry(
    client, query_vector, partition_numbers, top_k, collection_name="cids_fps", 
    job_id="default", max_retries=3, retry_delay=30
):
    """
    Search for similar vectors in specified partitions with retry logic.
    
    Args:
        client: MilvusClient instance
        query_vector: The query vector to search for
        partition_numbers: List of partition numbers to search in
        top_k: Number of results to return
        collection_name: Name of the Milvus collection
        job_id: Unique identifier for this job (for logging purposes)
        max_retries: Maximum number of retry attempts
        retry_delay: Delay between retry attempts in seconds
        
    Returns:
        List of search results
    """
    logger = get_progress_logger(job_id)
    partition_names = [f"cluster_{num}" for num in partition_numbers]
    
    for attempt in range(max_retries):
        try:
            logger.info(f"Searching in partitions: {', '.join(partition_names)} (attempt {attempt + 1}/{max_retries})")
            
            # Release collection first to ensure clean state
            try:
                client.release_collection(collection_name)
                logger.debug(f"Released collection: {collection_name}")
                # Small delay to allow for proper cleanup
                time.sleep(0.5)
            except Exception as e:
                logger.warning(f"Failed to release collection (continuing anyway): {str(e)}")

            # Load only the specific partitions we need
            logger.info(f"Loading partitions: {partition_names} for collection: {collection_name}")
            client.load_partitions(collection_name, partition_names=partition_names)
            logger.info(f"Partitions {partition_names} loaded successfully")
            
            # Small delay to ensure partitions are fully loaded
            time.sleep(1.0)
            
            # Perform search
            search_params = {
                "metric_type": "JACCARD",
                "params": {"nprobe": 10},
            }
            logger.info(f"Searching for top {top_k} results in collection: {collection_name}")
            results = client.search(
                collection_name=collection_name,
                data=[query_vector],
                anns_field="vector",
                search_params=search_params,
                limit=top_k,
                partition_names=partition_names,
            )
            logger.info(f"Search completed successfully in partitions: {partition_names}")
            # Immediately release partitions after successful search
            try:
                client.release_partitions(collection_name, partition_names=partition_names)
                logger.debug(f"Released partitions: {partition_names}")
            except Exception as e:
                logger.warning(f"Failed to release partitions (continuing anyway): {str(e)}")
            
            return results[0]
            
        except Exception as e:
            logger.error(f"Search attempt {attempt + 1} failed: {str(e)}")
            
            # Always try to clean up on error
            try:
                client.release_partitions(collection_name, partition_names=partition_names)
                logger.debug("Cleaned up partitions after error")
            except Exception as cleanup_error:
                logger.warning(f"Failed to cleanup partitions after error: {str(cleanup_error)}")
                
            # Also try to release the entire collection to ensure clean state
            try:
                client.release_collection(collection_name)
                logger.debug("Released entire collection for clean state")
            except Exception as release_error:
                logger.warning(f"Failed to release collection: {str(release_error)}")
            
            if attempt < max_retries - 1:
                logger.info(f"Waiting {retry_delay} seconds before retry...")
                time.sleep(retry_delay)
                
                # For connection-related errors, don't try to reconnect here since
                # the client is managed by the calling function
                if "connection" in str(e).lower() or "timeout" in str(e).lower():
                    logger.warning("Connection-related error detected - will be handled by caller")
            else:
                logger.error(f"All {max_retries} search attempts failed for partitions: {partition_names}")
                raise Exception(f"Failed to search partitions after {max_retries} attempts: {str(e)}")


def search_partitions(
    client, query_vector, partition_numbers, top_k, collection_name="cids_fps", job_id="default"
):
    """
    Search for similar vectors in specified partitions.
    
    Args:
        client: MilvusClient instance
        query_vector: The query vector to search for
        partition_numbers: List of partition numbers to search in
        top_k: Number of results to return
        collection_name: Name of the Milvus collection
        job_id: Unique identifier for this job (for logging purposes)
        
    Returns:
        List of search results
    """
    return search_partitions_with_retry(
        client, query_vector, partition_numbers, top_k, collection_name, job_id
    )


def run_search(query_vector, heapnum, job_id="default"):
    """
    Run a search against the Milvus database with retry logic.
    
    Args:
        query_vector: The query vector to search for
        heapnum: Number of results to return
        job_id: Unique identifier for this job (for logging purposes)
        
    Returns:
        Dictionary mapping CIDs to similarity scores
    """
    logger = get_progress_logger(job_id)
    top_k = heapnum
    candidate_list = []

    num_partitions = 120
    partitions_per_batch = int(os.getenv("PARTITION_DIVISION", 10))

    logger.info("Connecting to Milvus Client...")

    client = None
    try:
        client = connect_to_milvus_with_retry(job_id)
        logger.info("Beginning Vector Search.")

        total_batches = (num_partitions + partitions_per_batch - 1) // partitions_per_batch
        successful_batches = 0
        
        for i in range(0, num_partitions, partitions_per_batch):
            partition_batch = range(i, min(i + partitions_per_batch, num_partitions))
            batch_num = (i // partitions_per_batch) + 1
            
            try:
                logger.info(f"Processing batch {batch_num}/{total_batches}: partitions {list(partition_batch)}")
                results = search_partitions(client, query_vector, partition_batch, top_k=top_k, job_id=job_id)

                # For the search results of each partition, keep a list of the top_k best matches
                for hit in results:
                    id = str(hit["id"])
                    distance = 1 - hit["distance"]  # 0 is worst, 1 is best
                    if len(candidate_list) < top_k:
                        candidate_list.append((distance, id))
                    else:
                        # Find index of the worst match (largest distance, i.e., worst similarity)
                        worst_idx, worst_val = min(enumerate(candidate_list), key=lambda x: x[1][0])
                        if distance > worst_val[0]:
                            candidate_list[worst_idx] = (distance, id)

                successful_batches += 1
                percent_complete = (i + partitions_per_batch) / num_partitions * 100
                logger.info(f"Vector Search progress: {percent_complete:.2f}% complete ({successful_batches}/{total_batches} batches successful)")
                
            except Exception as e:
                logger.error(f"Failed to process batch {batch_num}/{total_batches}: {str(e)}")
                logger.warning(f"Continuing with remaining batches...")
                
                # If we have connection issues, try to reconnect for next batch
                if "connection" in str(e).lower() or "timeout" in str(e).lower():
                    logger.info("Connection issue detected, attempting to reconnect for next batch...")
                    try:
                        # Clean up old client first
                        if client:
                            cleanup_milvus_client(client, job_id)
                        client = connect_to_milvus_with_retry(job_id, max_retries=2, retry_delay=10)
                    except Exception as reconnect_error:
                        logger.error(f"Failed to reconnect: {str(reconnect_error)}")

        if successful_batches == 0:
            raise Exception("No batches were processed successfully")
        
        logger.info(f"Vector search completed. Successfully processed {successful_batches}/{total_batches} batches.")
        
        # Sort by distance (largest/best match first)
        candidate_list.sort(key=lambda x: x[0], reverse=True)
        return candidate_list
        
    finally:
        # Always cleanup the client connection
        if client:
            cleanup_milvus_client(client, job_id)


def runMilvus(query, heapnum, job_id="default"):
    """
    Main function to run a Milvus search using a query with retry logic.
    
    Args:
        query: Query string (SMILES or CID)
        heapnum: Number of results to return
        job_id: Unique identifier for this job (for logging purposes)
        
    Returns:
        Dictionary mapping CIDs to similarity scores
    """
    logger = get_progress_logger(job_id)
    # Ensure directories exist for this job
    ensure_directories_exist(job_id)
    max_retries = 3
    retry_delay = 30

    # Get the app state (singleton or passed in)
    try: 
        milvus_fp_cache = current_app.state.milvus_fp_cache
    except:
        milvus_fp_cache = None

    for attempt in range(max_retries):
        try:
            logger.info(f"Starting Milvus search with query: {query} (attempt {attempt + 1}/{max_retries})")
            # Convert the query FP into a Milvus searchable format
            fp_bytes = np.packbits(query).tobytes()
            cache_key = (fp_bytes, heapnum)
            if milvus_fp_cache and cache_key in milvus_fp_cache:
                logger.info("Returning cached Milvus search result.")
                logger.debug(f"Cache hit for key: {cache_key}. Cached result: {milvus_fp_cache[cache_key]}")
                milvus_fp_cache.move_to_end(cache_key)
                # Return a copy to prevent in-place mutation of the cached value
                return milvus_fp_cache[cache_key].copy()
            logger.debug(f"Cache miss for key: {cache_key}. Running new search.")
            results = run_search(fp_bytes, heapnum, job_id)
            logger.info("Milvus search completed successfully.")
            logger.debug(f"Search result for key {cache_key}: {results}")
            if results and len(results) > 0:
                # Only overwrite cache if new result is larger than existing
                if milvus_fp_cache and (cache_key not in milvus_fp_cache or len(results) > len(milvus_fp_cache[cache_key])):
                    logger.debug(f"Caching result for key {cache_key} (length: {len(results)}).")
                    milvus_fp_cache[cache_key] = results.copy()
                    if len(milvus_fp_cache) > 10:
                        oldest_key, oldest_val = milvus_fp_cache.popitem(last=False)
                        logger.debug(f"Evicted oldest cache entry: {oldest_key}")
                else:
                    if milvus_fp_cache: logger.debug(f"Not overwriting cache for key {cache_key} (existing length: {len(milvus_fp_cache[cache_key])}, new length: {len(results)}).")
            else:
                logger.debug(f"Not caching result for key {cache_key} (empty or None).")
            return results
        except Exception as e:
            logger.error(f"Milvus search attempt {attempt + 1} failed: {str(e)}")
            if attempt < max_retries - 1:
                logger.info(f"Waiting {retry_delay} seconds before retry...")
                time.sleep(retry_delay)
                # Force garbage collection to help with memory management
                import gc
                gc.collect()
                logger.debug("Forced garbage collection to free memory")
            else:
                logger.error(f"All {max_retries} attempts failed for Milvus search")
                raise Exception(f"Milvus search failed after {max_retries} attempts: {str(e)}")
    # This should never be reached, but just in case
    raise Exception("Unexpected error in runMilvus function")


def cleanup_milvus_client(client, job_id="default", collection_name="cids_fps"):
    """
    Properly cleanup Milvus client resources to prevent memory leaks.
    
    Args:
        client: MilvusClient instance to cleanup
        job_id: Unique identifier for this job (for logging purposes)
        collection_name: Name of the collection to release
    """
    logger = get_progress_logger(job_id)
    try:
        # Release all collections and partitions
        logger.info("Cleaning up Milvus client resources...")
        
        # Try to release the entire collection (this releases all partitions)
        try:
            client.release_collection(collection_name)
            logger.info(f"Released collection: {collection_name}")
        except Exception as e:
            logger.warning(f"Could not release collection {collection_name}: {str(e)}")
        
        # Close the client connection if it has a close method
        try:
            client.close()
        except Exception as e:
            logger.warning(f"Could not close client connection: {str(e)}")
            
        logger.info("Milvus client cleanup completed")
        
    except Exception as e:
        logger.error(f"Error during Milvus client cleanup: {str(e)}")
