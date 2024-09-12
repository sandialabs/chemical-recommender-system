# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import os
import logging
from App import create_app

# Configure Flask app logging
logging.basicConfig(
    filename="logs/app.log",
    filemode="w",  # Overwrite the log file each time
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

# Create a logger for this module
logger = logging.getLogger(__name__)

app = create_app()


def runFlask():
    """
    Starts the Flask server on the port specified in the environment variable CRS_PORT.
    If CRS_PORT is not set, it defaults to port 5005.
    """
    port = int(os.getenv("CRS_PORT", 5005))
    logger.info(f"Starting Flask server on port {port}")
    try:
        app.run(host="0.0.0.0", port=port)
    except Exception as e:
        logger.error(f"Error starting Flask server: {e}")
        raise


if __name__ == "__main__":
    runFlask()
