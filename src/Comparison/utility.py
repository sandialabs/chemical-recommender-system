# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import os
import sys
from pathlib import Path
from utils.logging_config import get_logger

logger = get_logger(__name__)

def ensure_directories_exist():
    """
    Create all required directories for the application to run properly.
    This prevents "Cannot save file into a non-existent directory" errors.
    """
    # Define all directories needed by the application
    directories = [
        "src/Comparison/LocalIO",
        "src/App/static/LocalIO",
        "src/App/static/jobs",
        "src/Comparison/Work",
        "src/Comparison/Output",
        "tmp"
    ]
    
    # Get the project root directory
    root_dir = Path(__file__).parent.parent.parent
    
    # Check and create each directory
    for directory in directories:
        try:
            dir_path = root_dir / directory
            dir_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Ensured directory exists: {dir_path}")
        except Exception as e:
            logger.error(f"Failed to create directory {directory}: {e}")
            
    # Also ensure job directory structure has proper permissions
    try:
        jobs_dir = root_dir / "src/App/static/jobs"
        jobs_dir.mkdir(parents=True, exist_ok=True)
        # Make sure directory has write permissions
        os.chmod(jobs_dir, 0o755)
    except Exception as e:
        logger.error(f"Failed to set permissions on jobs directory: {e}")
