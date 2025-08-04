# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np
import pubchempy as pcp


def get_compound_name(cid):
    """
    Get the compound name using PubChem CID.

    Inputs:
    - cid: The chemical identifier (CID).

    Returns:
    - name: The IUPAC name or a synonym of the compound, or "CID: <cid>" if not available.
    """
    if pcp is None:
        return "CID: " + str(cid)
        
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


def graph2D(qcid, querysmiles=None, logger=None):
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
    
    try:
        # Read the CSV file into a DataFrame
        df = pd.read_csv(csv_file)
        logger.debug(f"graph2D: Read CSV with shape {df.shape}")
        logger.debug(f"graph2D: CSV columns: {list(df.columns)}")
        
        if df.empty:
            logger.debug("graph2D: DataFrame is empty, cannot generate graph")
            return
            
        # The CSV structure from Controller.py is:
        # ['CID', 'Overall', 'Fingerprint', 'Molecular', 'Thermophysical', 'Toxicity', 'Synthetic']
        # Based on produceCSV: [CID, val[0], val[1], val[2], val[3], val[4], val[5]]
        # Where comp_dict[cid] = [overall, fingerprint/structural, molecular_weight, thermophysical, toxicity, SA_score]
        
        # Define the exact column names based on the actual CSV structure
        required_columns = {
            'cid_col': 'CID',
            'structural_col': 'Fingerprint',     # val[1] - Fingerprint/structural (Tanimoto) similarity
            'thermophysical_col': 'Thermophysical',  # val[3] - Thermophysical similarity  
            'toxicity_col': 'Toxicity',          # val[4] - Toxicity score
            'synthetic_col': 'Synthetic'         # val[5] - Synthetic availability (SA score)
        }
        
        # Verify all required columns exist
        missing_cols = []
        for col_name, col_key in required_columns.items():
            if col_key not in df.columns:
                missing_cols.append(col_key)
        
        if missing_cols:
            logger.debug(f"graph2D: Missing required columns: {missing_cols}")
            return
            
        logger.debug(f"graph2D: Using columns - Structural: {required_columns['structural_col']}, Thermophysical: {required_columns['thermophysical_col']}, Toxicity: {required_columns['toxicity_col']}, Synthetic: {required_columns['synthetic_col']}")
        
        # Convert DataFrame to records for processing
        records = df.to_dict(orient="records")
        logger.debug(f"graph2D: Processing {len(records)} data points")
        
        if len(records) < 1:
            logger.debug("graph2D: No data points available")
            return
        
    except Exception as e:
        logger.debug(f"graph2D: Error reading CSV: {e}")
        return

    # Use 'agg' backend for Matplotlib to avoid GUI issues
    matplotlib.use("agg")
    
    try:
        # Extract data points for plotting
        x_values = []
        y_values = []
        colors = []
        sizes = []
        
        # Process each record
        for i, record in enumerate(records):
            try:
                # Get structural similarity (X-axis)
                x = pd.to_numeric(record.get(required_columns['structural_col'], 0), errors='coerce')
                if pd.isna(x):
                    x = 0.0
                    
                # Get thermophysical similarity (Y-axis)
                y = pd.to_numeric(record.get(required_columns['thermophysical_col'], 0), errors='coerce')
                if pd.isna(y):
                    y = 0.0
                    
                # Get toxicity value (color)
                toxicity = pd.to_numeric(record.get(required_columns['toxicity_col'], 0.5), errors='coerce')
                if pd.isna(toxicity):
                    toxicity = 0.5
                    
                # Get synthetic availability (size)
                synthetic = pd.to_numeric(record.get(required_columns['synthetic_col'], 0.5), errors='coerce')
                if pd.isna(synthetic) or synthetic < 0:
                    synthetic = 0.5
                # Clamp synthetic values to reasonable range (0-10, SA scores are typically 1-10)
                synthetic = max(0.0, min(10.0, synthetic))
                
                x_values.append(float(x))
                y_values.append(float(y))
                colors.append(float(toxicity))
                
                # Calculate point size based on synthetic availability
                # SA scores range from 1-10, with higher = more synthesizable
                base_size = 100   # Base size for points
                max_size = 250   # Maximum size for highly synthesizable compounds
                size = base_size + (float(synthetic)**2) * (max_size - base_size)

                # Ensure size is always positive and finite
                if pd.isna(size) or size <= 0:
                    size = base_size
                sizes.append(size)
                
            except (ValueError, TypeError) as e:
                logger.debug(f"graph2D: Error processing record {i}: {e}")
                # Add default values to maintain array consistency
                x_values.append(0.0)
                y_values.append(0.0)
                colors.append(0.5)
                sizes.append(25.0)
        
        if not x_values:
            logger.debug("graph2D: No valid data points to plot")
            return
            
        logger.debug(f"graph2D: Successfully processed {len(x_values)} valid data points")
        
        # Create a custom colormap for toxicity values
        cmap = mcolors.LinearSegmentedColormap.from_list(
            "custom_cmap", [(0, "green"), (0.5, "yellow"), (1.0, "red")]
        )
        
        # Create the main plot
        fig, ax = plt.subplots(figsize=(10, 8))
        logger.debug(f"sizes: {sizes}, colors: {colors}, x_values: {x_values}, y_values: {y_values}")
        # Plot all points
        scatter = ax.scatter(x_values, y_values, c=colors, s=sizes, 
                           cmap=cmap, marker="o", edgecolors="black", alpha=0.7)
        
        # Highlight top matches if we have enough points
        if len(records) > 1:
            # Mark the top 3 matches (excluding query if it's the first one)
            start_idx = 1 if len(records) > 1 else 0
            for i, label in enumerate(["#1 Match", "#2 Match", "#3 Match"]):
                idx = start_idx + i
                if idx < len(records) and idx < len(x_values):
                    ax.scatter(x_values[idx], y_values[idx], c=colors[idx], s=sizes[idx],
                             cmap=cmap, marker="o", edgecolors="blue", linewidth=3,
                             label=label)
                    # Add annotation
                    ax.annotate(label, (x_values[idx], y_values[idx]), 
                              xytext=(5, 5), textcoords='offset points',
                              ha='left', va='bottom', fontsize=9,
                              bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
        
        # Set labels for the axes
        ax.set_xlabel("Structural Similarity")
        ax.set_ylabel("Thermophysical Similarity")
        
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
        
        # Add a colorbar to the plot
        colorbar = plt.colorbar(scatter, label="Predicted Toxicity", ax=ax)
        
        # Mark the toxicity of the query compound on the colorbar if available
        if colors:
            query_toxicity = colors[0] if colors else 0
            colorbar.ax.axhline(y=query_toxicity, color="blue", linestyle="-", linewidth=2)
        
        # Add explanatory text below the plot with proper spacing
        txt = "The radius of each point represents the Synthetic Availability of the candidate. The larger the point, the greater the availability of the candidate. Additionally, the blue line marked on the colorbar represents the toxicity of the query."
        plt.figtext(0.5, 0.01, txt, wrap=True, ha="center", va="bottom", fontsize=10)
        
        # Save the plot as a PNG file with adjusted layout
        plt.tight_layout()
        # Add bottom margin to prevent overlap between axis labels and caption
        plt.subplots_adjust(bottom=0.15)
        plt.savefig("src/App/static/LocalIO/graph.png", bbox_inches="tight", dpi=150)
        plt.close()
        logger.debug("graph2D: Successfully saved 2D graph")
        
    except Exception as e:
        logger.debug(f"graph2D: Error during graph generation: {e}")
        import traceback
        traceback.print_exc()
        # Create a simple fallback image
        try:
            plt.figure(figsize=(8, 6))
            plt.text(0.5, 0.5, f"Graph generation failed\nError: {str(e)}", 
                    ha='center', va='center', transform=plt.gca().transAxes, fontsize=12)
            plt.title("CRS Results - Graph Generation Error")
            plt.savefig("src/App/static/LocalIO/graph.png", bbox_inches="tight")
            plt.close()
            logger.debug("graph2D: Created fallback error graph")
        except Exception as fallback_error:
            logger.debug(f"graph2D: Could not create fallback graph: {fallback_error}")
        return


def graph3D(qcid, querysmiles=None, logger=None):
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
    
    try:
        # Read the CSV file into a DataFrame
        df = pd.read_csv(csv_file)
        logger.debug(f"graph3D: Read CSV with shape {df.shape}")
        
        if df.empty:
            logger.debug("graph3D: DataFrame is empty, cannot generate graph")
            return
            
        # The CSV structure from Controller.py is:
        # ['CID', 'Overall', 'Fingerprint', 'Molecular', 'Thermophysical', 'Toxicity', 'Synthetic']
        
        # Define the exact column names for 3D plotting
        required_columns = {
            'cid_col': 'CID',
            'structural_col': 'Fingerprint',     # X-axis: Fingerprint/structural (Tanimoto) similarity
            'thermophysical_col': 'Thermophysical',  # Y-axis: Thermophysical similarity
            'molecular_col': 'Molecular',        # Z-axis: Molecular weight similarity
            'toxicity_col': 'Toxicity',          # Color: Toxicity score
            'synthetic_col': 'Synthetic'         # Size: Synthetic availability (SA score)
        }
        
        # Verify all required columns exist
        missing_cols = []
        for col_name, col_key in required_columns.items():
            if col_key not in df.columns:
                missing_cols.append(col_key)
        
        if missing_cols:
            logger.debug(f"graph3D: Missing required columns: {missing_cols}")
            return
            
        logger.debug(f"graph3D: Using columns - X: {required_columns['structural_col']}, Y: {required_columns['thermophysical_col']}, Z: {required_columns['molecular_col']}, Color: {required_columns['toxicity_col']}, Size: {required_columns['synthetic_col']}")
        
        # Convert DataFrame to records for processing
        records = df.to_dict(orient="records")
        logger.debug(f"graph3D: Processing {len(records)} data points")
        
        if len(records) < 1:
            logger.debug("graph3D: No data points available")
            return
        
    except Exception as e:
        logger.debug(f"graph3D: Error reading CSV: {e}")
        return

    # Use 'agg' backend for Matplotlib to avoid GUI issues
    matplotlib.use("agg")
    
    try:
        # Extract data points for plotting
        x_values = []
        y_values = []
        z_values = []
        colors = []
        sizes = []
        
        # Process each record
        for i, record in enumerate(records):
            try:
                # Get structural similarity (X-axis)
                x = pd.to_numeric(record.get(required_columns['structural_col'], 0), errors='coerce')
                if pd.isna(x):
                    x = 0.0
                    
                # Get thermophysical similarity (Y-axis)
                y = pd.to_numeric(record.get(required_columns['thermophysical_col'], 0), errors='coerce')
                if pd.isna(y):
                    y = 0.0
                    
                # Get molecular weight similarity (Z-axis)
                z = pd.to_numeric(record.get(required_columns['molecular_col'], 0), errors='coerce')
                if pd.isna(z):
                    z = 0.0
                    
                # Get toxicity value (color)
                toxicity = pd.to_numeric(record.get(required_columns['toxicity_col'], 0.5), errors='coerce')
                if pd.isna(toxicity):
                    toxicity = 0.5
                
                # Get synthetic availability (size)
                synthetic = pd.to_numeric(record.get(required_columns['synthetic_col'], 5.0), errors='coerce')
                if pd.isna(synthetic) or synthetic < 0:
                    synthetic = 5.0
                # Clamp synthetic values to reasonable range (0-10, SA scores are typically 1-10)
                synthetic = max(0.0, min(10.0, synthetic))
                
                x_values.append(float(x))
                y_values.append(float(y))
                z_values.append(float(z))
                colors.append(float(toxicity))
                
                # Calculate size with proper scaling for 3D
                # SA scores range from 1-10, with higher = more synthesizable
                base_size = 50   # Base size for 3D points
                max_size = 200   # Maximum size for highly synthesizable compounds
                size = base_size + (float(synthetic)) * (max_size - base_size)
                
                # Ensure size is always positive and finite
                if pd.isna(size) or size <= 0:
                    size = base_size
                sizes.append(size)
                
            except (ValueError, TypeError) as e:
                logger.debug(f"graph3D: Error processing record {i}: {e}")
                # Add default values to maintain array consistency
                x_values.append(0.0)
                y_values.append(0.0)
                z_values.append(0.0)
                colors.append(0.5)
                sizes.append(75.0)
        
        if not x_values:
            logger.debug("graph3D: No valid data points to plot")
            return
            
        logger.debug(f"graph3D: Successfully processed {len(x_values)} valid data points")
        
        # Create a custom colormap for toxicity values
        cmap = mcolors.LinearSegmentedColormap.from_list(
            "custom_cmap", [(0, "green"), (0.5, "yellow"), (1.0, "red")]
        )
        
        # Create a 3D plot
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection="3d")
        
        # Plot all points
        scatter = ax.scatter(x_values, y_values, z_values, c=colors, s=sizes,
                           cmap=cmap, marker="o", edgecolors="black", alpha=0.7)
        
        # Highlight top matches if we have enough points
        if len(records) > 1:
            # Mark the top 3 matches (excluding query if it's the first one)
            start_idx = 1 if len(records) > 1 else 0
            for i, label in enumerate(["#1 Match", "#2 Match", "#3 Match"]):
                idx = start_idx + i
                if idx < len(records) and idx < len(x_values):
                    ax.scatter(x_values[idx], y_values[idx], z_values[idx], 
                             c=colors[idx], s=sizes[idx], cmap=cmap, marker="o", 
                             edgecolors="blue", linewidth=3, label=label)
        
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
        
        # Add explanatory text below the plot
        txt = "The radius of each point represents the Synthetic Availability of the candidate. The larger the point, the greater the availability of the candidate."
        plt.figtext(0.5, 0.02, txt, wrap=True, ha="center", va="center", fontsize=10)
        
        # Save the plot as a PNG file
        plt.tight_layout()
        plt.savefig("src/App/static/LocalIO/graph_3d.png", bbox_inches="tight", dpi=150)
        plt.close()
        logger.debug("graph3D: Successfully saved 3D graph")
        
    except Exception as e:
        logger.debug(f"graph3D: Error during 3D graph generation: {e}")
        import traceback
        traceback.print_exc()
        # Create a simple fallback image
        try:
            plt.figure(figsize=(8, 6))
            plt.text(0.5, 0.5, f"3D Graph generation failed\nError: {str(e)}", 
                    ha='center', va='center', transform=plt.gca().transAxes, fontsize=12)
            plt.title("CRS Results - 3D Graph Generation Error")
            plt.savefig("src/App/static/LocalIO/graph_3d.png", bbox_inches="tight")
            plt.close()
            logger.debug("graph3D: Created fallback error graph")
        except Exception as fallback_error:
            logger.debug(f"graph3D: Could not create fallback graph: {fallback_error}")
        return


def graphResults(qcid, querysmiles=None, logger=None):
    """
    Generate both 2D and 3D graphs of the results.

    Inputs:
    - qcid: The query chemical identifier (CID).
    - querysmiles: The query SMILES string (optional).

    Returns:
    - None (saves the plots as PNG files).
    """
    # Generate 2D graph
    graph2D(qcid, querysmiles, logger)
    # Generate 3D graph
    graph3D(qcid, querysmiles, logger)
