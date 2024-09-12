# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

from flask import Blueprint, render_template

contact_bp = Blueprint("contact", __name__)


@contact_bp.route("/contact", methods=["GET"])
def contact():
    """
    Renders the contact page.

    This function handles GET requests to the /contact route and renders the contact.html template.

    Inputs:
    - None (uses no input variables directly)

    Returns:
    - Renders the 'contact.html' template.
    """
    return render_template("contact.html")
