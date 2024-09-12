# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

from flask import Blueprint, render_template

main_bp = Blueprint("main", __name__)


@main_bp.route("/", methods=["GET", "POST"])
def index():
    """
    Renders the main index page.

    This function handles both GET and POST requests to the root route ("/") and renders the index.html template.

    Inputs:
    - None (uses no input variables directly)

    Returns:
    - Renders the 'index.html' template.
    """
    return render_template("index.html")
