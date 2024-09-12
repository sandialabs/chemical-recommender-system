# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import pandas as pd
from pymilvus import (
    connections,
    Collection,
    FieldSchema,
    CollectionSchema,
    DataType,
    utility,
)
from tqdm.auto import tqdm
import zipfile
import os


def convert_bit_string_to_bytes(bit_string):
    try:
        num_bytes = (len(bit_string) + 7) // 8
        return int(bit_string, 2).to_bytes(num_bytes, byteorder="big")
    except ValueError:
        if bit_string != "nan":
            print(f"Error converting bit_string: {bit_string}")
        return None


def process_and_insert_chunk(collection, partition_name, chunk):
    chunk["fps"] = chunk["fps"].astype(str)
    chunk["fps"] = [convert_bit_string_to_bytes(x) for x in chunk["fps"]]
    df_filtered = chunk.dropna(subset=["fps"])
    cids = df_filtered["cid"].tolist()
    vectors = df_filtered["fps"].tolist()
    collection.insert([cids, vectors], partition_name=partition_name)


def read_csv_and_insert(zip_file_path, collection, partition_name, chunk_size=10000):
    with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
        csv_file_name = zip_ref.namelist()[0]  # Assuming each zip contains a single CSV
        with zip_ref.open(csv_file_name) as myzip:
            reader = pd.read_csv(myzip, chunksize=chunk_size)
            for chunk in tqdm(
                reader, desc=f"Processing and Inserting into {partition_name}"
            ):
                process_and_insert_chunk(collection, partition_name, chunk)


def main():
    # Connect to Milvus
    print("Connecting to Milvus...")
    connections.connect(host="localhost", port="19530")
    collection_name = "cids_fps"  # USE HERE TO RENAME

    # Drop the collection if it exists
    if utility.has_collection(collection_name):
        collection = Collection(name=collection_name)
        collection.drop()
        print(f"Collection '{collection_name}' has been dropped.")

    # Define schema
    id_field = FieldSchema(name="cid", dtype=DataType.INT64, is_primary=True)
    vector_field = FieldSchema(name="vector", dtype=DataType.BINARY_VECTOR, dim=2048)
    schema = CollectionSchema(
        fields=[id_field, vector_field], description="Chemical Compounds"
    )

    # Create collection in Milvus
    print("Creating collection in Milvus...")
    collection = Collection(name=collection_name, schema=schema)

    # Directory containing cluster zips
    clusters_dir = "models/clusters/"

    # Process each cluster file
    for i in range(120):  # Assuming cluster IDs range from 0 to 19
        cluster_zip_path = os.path.join(clusters_dir, f"cluster{i}.zip")
        partition_name = f"cluster_{i}"
        collection.create_partition(partition_name=partition_name)
        read_csv_and_insert(cluster_zip_path, collection, partition_name)

    # Creating index
    print("Creating index...")
    collection.create_index(
        field_name="vector",
        index_params={
            "index_type": "BIN_IVF_FLAT",
            "metric_type": "JACCARD",
            "params": {"nlist": 1024},
        },
    )

    print("Process completed.")


if __name__ == "__main__":
    main()
