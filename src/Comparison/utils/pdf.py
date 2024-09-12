# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause
import csv
import pubchempy as pcp
from reportlab.lib.pagesizes import letter, landscape
from reportlab.pdfgen import canvas
from io import BytesIO
from fpdf import FPDF
from pypdf import PdfWriter


def fetchPubchemInfo(cids):
    """
    Fetch compound properties from PubChem for a list of CIDs.

    Inputs:
    - cids: A list of chemical identifiers (CIDs).

    Returns:
    - properties: A list of dictionaries containing compound properties.
    """
    cidarr = [int(cid) for cid in cids if cid != "-1"]
    try:
        properties = pcp.get_properties(["CanonicalSMILES", "IUPACName"], cidarr)

        for prop in properties:
            if not prop.get("IUPACName"):
                synonyms = pcp.get_synonyms(prop["CID"])
                if synonyms:
                    prop["IUPACName"] = synonyms[0]["Synonym"][0]
                else:
                    prop["IUPACName"] = "cid:" + str(prop["CID"])

        return properties
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def get_compound_name(properties, cid):
    """
    Get the compound name from the fetched properties.

    Inputs:
    - properties: A list of dictionaries containing compound properties.
    - cid: The chemical identifier (CID).

    Returns:
    - name: The IUPAC name or a synonym of the compound, or "CID: <cid>" if not available.
    """
    try:
        # Find the property dictionary for the given CID
        property_dict = next(item for item in properties if item["CID"] == cid)
        # Extract the IUPAC name from the properties
        name = property_dict.get("IUPACName", None)
        if name is None:
            # If IUPAC name is not available, fetch synonyms
            synonyms = property_dict.get("Synonym", [])
            # Use the first synonym if available, otherwise use CID
            name = synonyms[0] if synonyms else "CID: " + str(cid)
        return name
    except:
        # Return CID if any error occurs
        return "CID: " + str(cid)


def draw_string_with_dynamic_font(
    c, x, y, text, max_width, centered=True, initial_font_size=12, font_name="Helvetica"
):
    """
    Draw a centered string with a dynamically adjusted font size.

    Parameters:
    - c (Canvas): The ReportLab canvas object.
    - x (float): The x-coordinate for the center of the text.
    - y (float): The y-coordinate for the center of the text.
    - text (str): The text to be drawn.
    - max_width (float): The maximum allowed width for the text.
    - initial_font_size (int, optional): The initial font size. Default is 12.
    - font_name (str, optional): The font name. Default is "Helvetica".

    Returns:
    - None
    """
    # Set the initial font size
    font_size = initial_font_size
    c.setFont(font_name, font_size)

    # Calculate the width of the text
    text_width = c.stringWidth(text, font_name, font_size)

    # Reduce the font size until the text fits within the max_width
    while text_width > max_width and font_size > 1:
        font_size -= 1
        c.setFont(font_name, font_size)
        text_width = c.stringWidth(text, font_name, font_size)

    # Draw the text centered at (x, y)
    if centered:
        c.drawCentredString(x, y, text)
    else:
        c.drawString(x, y, text)


def createPDFImages(pubchem_cids, output_file="src/Comparison/LocalIO/cidinfo.pdf"):
    """
    Create a PDF with images of compounds from PubChem.

    Inputs:
    - pubchem_cids: A list of chemical identifiers (CIDs).
    - output_file: The output file path for the PDF (default is "src/Comparison/LocalIO/cidinfo.pdf").

    Returns:
    - None (saves the PDF to the specified output file).
    """
    # Create a buffer to hold the PDF data
    buffer = BytesIO()
    # Create a PDF canvas with landscape orientation
    pdf = canvas.Canvas(buffer, pagesize=landscape(letter))

    # Define positions for images on the PDF page
    bases = [
        [37.5, 550],
        [287.5, 550],
        [537.5, 550],
        [37.5, 280],
        [287.5, 280],
        [537.5, 280],
    ]
    num = 1
    first = True
    totalnum = 1

    # Fetch properties for all CIDs at once
    compound_properties = fetchPubchemInfo(pubchem_cids)

    # Iterate over the list of CIDs
    for cid in pubchem_cids:
        if cid == "-1" or cid == -1:
            # Reset counters if CID is -1
            totalnum = 1
            first = False
        else:
            if compound_properties:
                if num == 7:
                    # Start a new page if more than 6 images
                    num = 1
                    pdf.translate(0, -800)
                    pdf.showPage()

                # Get the compound name for the current CID
                name = get_compound_name(compound_properties, cid)
                if len(name) > 38:
                    name = name[:34] + "..."
                cid_val = cid
                # URL to fetch the compound image from PubChem
                image_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_val}/PNG"
                base = bases[num - 1]
                if first:
                    totalnum = "Query"
                # Draw the compound name and CID on the PDF
                pdf.drawString(base[0], base[1], f"{totalnum}.{name}")
                if first:
                    totalnum = 0
                pdf.drawString(base[0], base[1] - 20, f"CID: {cid_val}")
                # Draw the compound image on the PDF
                pdf.drawImage(image_url, base[0], base[1] - 250, width=220, height=220)

                num += 1
                totalnum += 1
            first = False

    # Save the PDF to the buffer
    pdf.save()

    # Write the buffer content to the output file
    with open(output_file, "wb") as f:
        f.write(buffer.getvalue())


def convert(
    source: str,
    destination: str,
    orientation="P",
    delimiter=",",
    align="C",
    size=8,
    headersize=9,
) -> None:
    """
    Convert a CSV file to a PDF table.

    Inputs:
    - source: The source CSV file path.
    - destination: The destination PDF file path.
    - orientation: The page orientation ("P" for portrait, "L" for landscape, default is "P").
    - delimiter: The delimiter used in the CSV file (default is ",").
    - align: The alignment for the table cells (default is "C" for center).
    - size: The font size for the data rows (default is 8).
    - headersize: The font size for the header row (default is 9).

    Returns:
    - None (saves the PDF to the specified destination).
    """
    # Validate font size parameters
    if not (isinstance(size, int) and isinstance(headersize, int)):
        raise Exception("Type Error: Font Size should be of int data type")

    # Create a PDF object with the specified orientation and page size
    PDF = FPDF(orientation, format="letter")
    PDF.add_page()

    # Read the CSV file
    with open(source) as CSV:
        data = [row for row in csv.reader(CSV, delimiter=delimiter)]
        header = data[0]  # Extract header row
        rows = data[1:]  # Extract data rows

    # Determine the maximum number of columns
    max_cols = len(header)
    for row in rows:
        if len(row) > max_cols:
            max_cols = len(row)

    # Extend header and rows to match the maximum number of columns
    header.extend([" "] * (max_cols - len(header)))
    for row in rows:
        row.extend([" "] * (max_cols - len(row)))

    # Set font for the header
    PDF.set_font("Helvetica", "B", size=headersize)

    # Calculate line heights and column width
    header_line_height = PDF.font_size * 4.15
    row_line_height = PDF.font_size * 5.25
    col_width = PDF.epw / max_cols

    # Function to add a row to the PDF table
    def add_table_row(cells, line_height, is_header=False):
        for cell in cells:
            x_before = PDF.get_x()
            y_before = PDF.get_y()
            cell_height = line_height

            # Calculate vertical offset for centering text
            if is_header:
                PDF.set_font("Helvetica", "B", size=headersize)
                text_height = PDF.get_string_width(cell) / col_width * headersize
            else:
                PDF.set_font("Helvetica", size=size)
                text_height = PDF.get_string_width(cell) / col_width * size

            vertical_offset = (cell_height - text_height) / 2

            # Draw the cell with vertical centering
            PDF.multi_cell(
                col_width,
                line_height,
                cell,
                align=align,
                border=1,
                max_line_height=PDF.font_size,
            )
            PDF.set_xy(x_before + col_width, y_before)
        PDF.set_y(y_before + line_height)

    # Add the header row to the PDF table
    add_table_row(header, header_line_height, is_header=True)

    # Set font for the data rows
    PDF.set_font("Helvetica", size=size)

    # Add data rows to the PDF table
    for row in rows:
        if PDF.get_y() + row_line_height > PDF.page_break_trigger:
            PDF.add_page()
            add_table_row(header, header_line_height, is_header=True)
        add_table_row(row, row_line_height)

    # Output the PDF to the specified destination
    PDF.output(destination)


def createPDFFirst(qval, qcid, params, weights, subfailed, containers):
    """
    Create the first page of the PDF report.

    Inputs:
    - qval: The query value.
    - qcid: The query chemical identifier (CID).
    - params: A list of parameters for the query.
    - weights: A list of weights for the query.
    - subfailed: A flag indicating if substructure searching failed.
    - containers: A list of container names for additional columns.

    Returns:
    - None (saves the PDF to "src/Comparison/LocalIO/first.pdf").
    """
    tarray = params[2]
    endpoints = ""

    # Extract thermophysical properties from parameters
    MP = tarray[0]
    if MP:
        endpoints += "MP, "
    BP = tarray[1]
    if BP:
        endpoints += "BP, "
    logP = tarray[2]
    if logP:
        endpoints += "logP, "
    Hlaw = tarray[3]
    if Hlaw:
        endpoints += "Hlaw, "
    VP = tarray[4]
    if VP:
        endpoints += "VP, "

    # Remove trailing comma and space
    endpoints = endpoints[:-2]

    # Check if all elements included
    incEle = params[3]
    if incEle == [""]:
        incEle = "None"
    if incEle == True:
        incEleStr = "Include all elements is on"
    else:
        incEleStr = "Including elements: " + str(incEle)

    # Extract SMARTS parameters
    smarts = params[4]
    smarts_num = params[5]
    if not smarts:
        smarts_string = "Substructure Searching is: off"
    else:
        if smarts_num:
            if not subfailed:
                smarts_string = (
                    "Substructure searching for "
                    + smarts
                    + " for "
                    + str(smarts_num)
                    + " appearance(s)"
                )
            else:
                smarts_string = (
                    "FAILED: Substructure searching for "
                    + smarts
                    + " for "
                    + str(smarts_num)
                    + " appearance(s)"
                )
        else:
            if not subfailed:
                smarts_string = (
                    "Substructure searching for " + smarts + " for >= 1 appearance(s)"
                )
            else:
                smarts_string = (
                    "FAILED: Substructure searching for "
                    + smarts
                    + " for >= 1 appearance(s)"
                )

    # Output file path for the first PDF page
    output_file = "src/Comparison/LocalIO/first.pdf"
    # Create a buffer to hold the PDF data
    buffer = BytesIO()
    # Create a PDF canvas with landscape orientation
    pdf = canvas.Canvas(buffer, pagesize=landscape(letter))
    pwidth, pheight = pdf._pagesize

    # Set font and draw the title
    draw_string_with_dynamic_font(
        pdf,
        pwidth / 2,
        pheight - 50,
        "CRS Result for: " + str(qval),
        550,
        True,
        18,
        "Helvetica-Bold",
    )
    pdf.line(pwidth / 2 - 300, pheight - 60, pwidth / 2 + 300, pheight - 60)

    # Draw settings and parameters on the PDF
    pdf.setFont("Helvetica-Bold", 14)
    pdf.drawString(85, pheight - 100, "Settings")
    pdf.setFont("Helvetica", 12)
    pdf.drawString(85, pheight - 120, incEleStr)
    pdf.drawString(85, pheight - 140, smarts_string)
    pdf.drawString(85, pheight - 160, f"Thermophysical Properties: {endpoints}")
    pdf.drawString(85, pheight - 180, f"Weightages: {weights}")
    pdf.drawString(85, pheight - 200, f"Model Images: {containers}")

    # Fetch and draw the query compound information if available
    if qcid != -1:
        compound_info = fetchPubchemInfo([qcid])
        if compound_info:
            name = get_compound_name(compound_info, qcid)
            cid_val = qcid
            image_url = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_val}/PNG"
            )
            pdf.setFont("Helvetica-Bold", 14)
            draw_string_with_dynamic_font(
                pdf,
                85,
                pheight - 240,
                f"QUERY: {name}",
                268,
                False,
                14,
                "Helvetica-Bold",
            )
            pdf.line(85, pheight - 245, 85 + 275, pheight - 245)
            pdf.drawImage(image_url, 85, pheight - 540, width=275, height=275)
    else:
        draw_string_with_dynamic_font(
            pdf,
            85,
            pheight - 240,
            "Query not found in PubChem (No Visualization Below)",
            320,
            False,
            14,
            "Helvetica-Bold",
        )
        pdf.line(85, pheight - 245, 430, pheight - 245)

    # Draw the 2D and 3D graphs on the PDF
    pdf.drawImage(
        "src/App/static/LocalIO/graph.png",
        pwidth - 370,
        pheight - 305,
        width=300,
        height=235,
    )
    pdf.drawImage(
        "src/App/static/LocalIO/graph_3d.png",
        pwidth - 370,
        pheight - 555,
        width=300,
        height=230,
    )

    # Save the PDF to the buffer
    pdf.save()

    # Write the buffer content to the output file
    with open(output_file, "wb") as f:
        f.write(buffer.getvalue())


def createPDFCSV():
    """
    Create a PDF from a CSV file.

    Inputs:
    - None

    Returns:
    - None (saves the PDF to "src/Comparison/LocalIO/csvinfo.pdf").
    """
    convert(
        "src/App/static/LocalIO/results.csv",
        "src/Comparison/LocalIO/csvinfo.pdf",
        orientation="L",
        size=9,
    )


def combinePDFs(
    pubchem_cids, queryval, qcid, params, weights, include, subfailed, containers=None
):
    """
    Combine multiple PDFs into a single PDF report.

    Inputs:
    - pubchem_cids: A list of chemical identifiers (CIDs).
    - queryval: The query value.
    - qcid: The query chemical identifier (CID).
    - params: A list of parameters for the query.
    - weights: A list of weights for the query.
    - include: A flag indicating whether to include additional PDFs.
    - subfailed: A flag indicating if substructure searching failed.
    - containers: A list of container names for additional columns (default is None).

    Returns:
    - None (saves the merged PDF to "src/App/static/LocalIO/report.pdf").
    """
    # Create individual PDFs
    createPDFImages(pubchem_cids)
    createPDFCSV()
    createPDFFirst(queryval, qcid, params, weights, subfailed, containers)

    # Initialize a PDF merger
    merger = PdfWriter()

    # List of PDFs to be merged
    if include:
        pdfs = [
            "src/App/static/assets/CRS1pagesum.pdf",
            "src/Comparison/LocalIO/first.pdf",
            "src/Comparison/LocalIO/csvinfo.pdf",
            "src/Comparison/LocalIO/cidinfo.pdf",
        ]
    else:
        pdfs = [
            "src/Comparison/LocalIO/first.pdf",
            "src/Comparison/LocalIO/csvinfo.pdf",
            "src/Comparison/LocalIO/cidinfo.pdf",
        ]

    # Append each PDF to the merger
    for pdf in pdfs:
        merger.append(pdf)

    # Write the merged PDF to the output file
    merger.write("src/App/static/LocalIO/report.pdf")
    merger.close()
