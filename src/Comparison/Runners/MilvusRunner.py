# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import time
import heapq
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pymilvus import MilvusClient
import logging

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


def smiles_to_fingerprint(smiles, radius=2, nbits=2048, useFeatures=True):
    """
    Generate a Morgan fingerprint from a SMILES string.

    Inputs:
    - smiles: A string representing the SMILES notation of a molecule.
    - radius: The radius of the Morgan fingerprint (default is 2).
    - nbits: The size of the fingerprint in bits (default is 2048).
    - useFeatures: Whether to use feature invariants (default is True).

    Returns:
    - arr: A NumPy array representing the Morgan fingerprint.
    """
    logger.debug(f"Generating fingerprint for SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius, nBits=nbits, useFeatures=useFeatures
    )
    # Convert the fingerprint to a NumPy array
    arr = np.zeros((1,))
    AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
    logger.debug(f"Generated fingerprint: {arr}")
    return arr


def search_partitions(
    client, query_vector, partition_numbers, top_k, collection_name="cids_fps"
):
    """
    Search for the most similar compounds to the query vector in specific partitions of the Milvus collection.

    Inputs:
    - client: An instance of the MilvusClient.
    - query_vector: The query vector to search for.
    - partition_numbers: A list of partition numbers to search within.
    - top_k: The number of top results to return.
    - collection_name: The name of the Milvus collection (default is "cids_fps").

    Returns:
    - results[0]: The search results for the query vector.
    """
    partition_names = [f"cluster_{num}" for num in partition_numbers]
    client.release_collection(collection_name)

    try:
        # Perform the search within the specified partitions
        client.load_partitions(collection_name, partition_names=partition_names)

        search_params = {
            "metric_type": "JACCARD",
            "params": {"nprobe": 10},
        }
        results = client.search(
            collection_name=collection_name,
            data=[query_vector],
            anns_field="vector",
            search_params=search_params,
            limit=top_k,
            partition_names=partition_names,  # Specify the partitions to search
        )
    finally:
        # Ensure partitions are released afterward or RAM will get filled and Docker will crash
        client.release_partitions(collection_name, partition_names=partition_names)

    return results[0]


def run_search(query_vector, heapnum):
    """
    Run a vector search across multiple partitions and return the top results.

    Inputs:
    - query_vector: The query vector to search for.
    - heapnum: The number of top results to return.

    Returns:
    - max_heap: A list of tuples containing the closest CIDs and their hit distances.
    """
    top_k = heapnum
    max_heap = []

    start_time = time.time()

    num_partitions = 120
    partitions_per_batch = 24

    client = MilvusClient("tcp://standalone:19530")

    logger.info("Successfully connected to Milvus Client. Beginning Vector Search.")

    total_batches = (num_partitions + partitions_per_batch - 1) // partitions_per_batch
    for i in range(0, num_partitions, partitions_per_batch):
        partition_batch = range(i, min(i + partitions_per_batch, num_partitions))
        results = search_partitions(client, query_vector, partition_batch, top_k=top_k)

        # For the search results of each partition, keep a max heap with the highest similarity values
        for hit in results:
            id = str(hit["id"])
            distance = 1 - hit["distance"]
            if len(max_heap) > top_k:
                heapq.heapreplace(max_heap, (distance, id))
            else:
                heapq.heappush(max_heap, (distance, id))

        percent_complete = (i + partitions_per_batch) / num_partitions * 100
        logger.info(f"Vector Search progress: {percent_complete:.2f}% complete")

    end_time = time.time()
    logger.info(
        f"Total time taken for vector search: {end_time - start_time:.2f} seconds"
    )

    # Return max_heap with closest CIDs from vector search and their hit distances
    return max_heap


def runMilvus(query, heapnum):
    """
    Convert the query fingerprint into a Milvus searchable format and run the search.

    Inputs:
    - query: The query fingerprint.
    - heapnum: The number of top results to return.

    Returns:
    - The result of the run_search function.
    """
    logger.info("============ Entering Milvus search ============")
    # Convert the query FP into a Milvus searchable format
    bit_string = query.ToBitString()
    fp_bytes = int(bit_string, 2).to_bytes((len(bit_string) + 7) // 8, byteorder="big")
    results = run_search(fp_bytes, heapnum)
    return results
