# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pubchempy as pcp
import pandas as pd
import numpy as np


def get_compound_name(cid):
    """
    Get the compound name using PubChem CID.

    Inputs:
    - cid: The chemical identifier (CID).

    Returns:
    - name: The IUPAC name or a synonym of the compound, or "CID: <cid>" if not available.
    """
    try:
        # Fetch properties of the compound using PubChem CID
        properties = pcp.get_properties(["IUPACName"], cid)[0]
        # Extract the IUPAC name from the properties
        name = properties.get("IUPACName", None)
        if name is None:
            # If IUPAC name is not available, fetch synonyms
            synonyms = pcp.get_synonyms(cid)[0].get("Synonym", [])
            # Use the first synonym if available, otherwise use CID
            name = synonyms[0] if synonyms else "CID: " + str(cid)
        return name
    except:
        # Return CID if any error occurs
        return "CID: " + str(cid)


def graph2D(qcid, querysmiles=None):
    """
    Plot a 2D graph of the results.

    Inputs:
    - qcid: The query chemical identifier (CID).
    - querysmiles: The query SMILES string (optional).

    Returns:
    - None (saves the plot as a PNG file).
    """
    # Path to the CSV file containing data
    csv_file = r"src/Comparison/LocalIO/data.csv"
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file)
    # Convert DataFrame to dictionary with lists as values
    dictionary = df.to_dict(orient="list")

    # Use 'agg' backend for Matplotlib to avoid GUI issues
    matplotlib.use("agg")
    # Extract keys from the dictionary
    keys = list(dictionary.keys())

    for key, val in dictionary.items():
        x = val[1]  # Structural similarity
        y = val[3]  # Thermophysical similarity
        size = val[5]  # Synthetic availability

        s = size * 3**8 / 20

        plt.scatter(x, y, s=size * 3**8 / 20, marker="o", edgecolors="black")

    ax = plt.gca()
    tick = np.diff(ax.get_yticks())[0]
    plt.clf()

    # First key corresponds to the query compound
    first_key = keys[0]
    # Extract toxicity value of the query compound
    qtox = dictionary[first_key][4]
    # Extract keys for the top matches
    second_key = keys[1] if len(keys) > 1 else None
    third_key = keys[2] if len(keys) > 2 else None
    fourth_key = keys[3] if len(keys) > 3 else None

    # Create a custom colormap for toxicity values
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "custom_cmap", [(0, "green"), (0.5, "yellow"), (1.0, "red")]
    )

    # Iterate over the dictionary items to plot each point
    for key, val in dictionary.items():
        x = val[1]  # Structural similarity
        y = val[3]  # Thermophysical similarity
        colorval = val[4]  # Toxicity value
        size = val[5]  # Synthetic availability
        color = cmap(colorval)  # Color based on toxicity

        s = size * 3**8 / 20
        # Ignore the first key for the query compound
        if key == first_key:
            continue
        elif key == second_key:
            # Plot the first match with specific marker and annotation
            plt.scatter(
                x,
                y,
                color=color,
                marker="o",
                s=s,
                label="#1 Match",
                edgecolors="black",
            )
            plt.annotate(
                "#1 Match",
                (x, y + tick * 0.3),
                annotation_clip=False,
                ha="center",
            )
        elif key == third_key:
            # Plot the second match with specific marker and annotation
            plt.scatter(
                x,
                y,
                color=color,
                marker="o",
                s=s,
                label="#2 Match",
                edgecolors="black",
            )
            plt.annotate(
                "#2 Match",
                (x, y + tick * 0.3),
                annotation_clip=False,
                ha="center",
            )
        elif key == fourth_key:
            # Plot the third match with specific marker and annotation
            plt.scatter(
                x,
                y,
                color=color,
                marker="o",
                s=s,
                label="#3 Match",
                edgecolors="black",
            )
            plt.annotate(
                "#3 Match",
                (x, y + tick * 0.3),
                annotation_clip=False,
                ha="center",
            )
        else:
            # Plot other points without specific labels
            plt.scatter(
                x, y, color=color, s=size * 3**8 / 20, marker="o", edgecolors="black"
            )

    # Set labels for the axes
    plt.xlabel("Structural Similarity")
    plt.ylabel("Thermophysical Similarity")

    # Set the title based on whether SMILES or CID is provided
    if querysmiles is None:
        gname = get_compound_name(qcid)
        if len(gname) > 20:
            gname = gname[:20] + "..."
        plt.title("CRS Results for: CID " + str(qcid) + ", " + gname)
    else:
        gname = querysmiles
        if len(gname) > 20:
            gname = gname[:20] + "..."
        plt.title("CRS Results for SMILES: " + gname)
    plt.legend()

    # Add explanatory text below the plot
    txt = "The radius of each point represents the Synthetic Availability of the candidate. The larger the point, the greater the availability of the candidate. Additionally, the blue line marked on the colorbar represents the toxicity of the query."
    plt.figtext(0.5, -0.05, txt, wrap=True, ha="center", va="center", fontsize=10)

    # Create a ScalarMappable for the colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array([])

    # Add a colorbar to the plot
    ax = plt.gca()
    cax = ax.figure.add_axes([0.92, 0.1, 0.03, 0.8])
    colorbar = plt.colorbar(
        sm, label="Predicted Toxicity", cax=cax, orientation="vertical"
    )
    # Mark the toxicity of the query compound on the colorbar
    colorbar.ax.axhline(y=qtox, color="blue", linestyle="-")

    # Save the plot as a PNG file
    plt.savefig("src/App/static/LocalIO/graph.png", bbox_inches="tight")
    plt.clf()


def graph3D(qcid, querysmiles=None):
    """
    Plot a 3D graph of the results.

    Inputs:
    - qcid: The query chemical identifier (CID).
    - querysmiles: The query SMILES string (optional).

    Returns:
    - None (saves the plot as a PNG file).
    """
    # Path to the CSV file containing data
    csv_file = r"src/Comparison/LocalIO/data.csv"
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file)
    # Convert DataFrame to dictionary with lists as values
    dictionary = df.to_dict(orient="list")

    # Use 'agg' backend for Matplotlib to avoid GUI issues
    matplotlib.use("agg")
    # Extract keys from the dictionary
    keys = list(dictionary.keys())
    # First key corresponds to the query compound
    first_key = keys[0]
    # Extract keys for the top matches
    second_key = keys[1] if len(keys) > 1 else None
    third_key = keys[2] if len(keys) > 2 else None
    fourth_key = keys[3] if len(keys) > 3 else None

    # Create a custom colormap for toxicity values
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "custom_cmap", [(0, "green"), (0.5, "yellow"), (1.0, "red")]
    )

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Iterate over the dictionary items to plot each point
    for key, val in dictionary.items():
        x = val[1]  # Structural similarity
        y = val[3]  # Thermophysical similarity
        z = val[2]  # Molecular weight similarity
        colorval = val[4]  # Toxicity value
        size = val[5]  # Synthetic availability
        color = cmap(colorval)  # Color based on toxicity

        # Ignore the first key for the query compound
        if key == first_key:
            continue
        elif key == second_key:
            # Plot the first match with specific marker
            ax.scatter(
                x,
                y,
                z,
                color=color,
                marker="o",
                s=size * 3**8 / 20,
                label="#1 Match",
                edgecolors="black",
            )
        elif key == third_key:
            # Plot the second match with specific marker
            ax.scatter(
                x,
                y,
                z,
                color=color,
                marker="o",
                s=size * 3**8 / 20,
                label="#2 Match",
                edgecolors="black",
            )
        elif key == fourth_key:
            # Plot the third match with specific marker
            ax.scatter(
                x,
                y,
                z,
                color=color,
                marker="o",
                s=size * 3**8 / 20,
                label="#3 Match",
                edgecolors="black",
            )
        else:
            # Plot other points without specific labels
            ax.scatter(
                x, y, z, color=color, s=size * 3**8 / 20, marker="o", edgecolors="black"
            )

    # Set labels for the axes
    ax.set_xlabel("Structural Similarity")
    ax.set_ylabel("Thermophysical Similarity")
    ax.set_zlabel("Molecular Weight Similarity")

    # Set the title based on whether SMILES or CID is provided
    if querysmiles is None:
        gname = get_compound_name(qcid)
        if len(gname) > 20:
            gname = gname[:20] + "..."
        ax.set_title("CRS Results for: CID " + str(qcid) + ", " + gname)
    else:
        gname = querysmiles
        if len(gname) > 20:
            gname = gname[:20] + "..."
        ax.set_title("CRS Results for SMILES: " + gname)
    ax.legend()

    # Add explanatory text below the plot
    txt = "The radius of each point represents the Synthetic Availability of the candidate. The larger the point, the greater the availability of the candidate."
    plt.figtext(0.5, -0.05, txt, wrap=True, ha="center", va="center", fontsize=10)

    # Save the plot as a PNG file
    plt.savefig("src/App/static/LocalIO/graph_3d.png", bbox_inches="tight")
    plt.clf()


def graphResults(qcid, querysmiles=None):
    """
    Generate both 2D and 3D graphs of the results.

    Inputs:
    - qcid: The query chemical identifier (CID).
    - querysmiles: The query SMILES string (optional).

    Returns:
    - None (saves the plots as PNG files).
    """
    # Generate 2D graph
    graph2D(qcid, querysmiles)
    # Generate 3D graph
    graph3D(qcid, querysmiles)
