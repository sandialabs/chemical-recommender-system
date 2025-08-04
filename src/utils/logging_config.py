# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import logging
import os

def setup_logging(reset=False):
    """
    Set up centralized logging configuration for the entire application.
    
    Args:
        reset (bool): If True, reset all existing loggers (for new runs)
                      If False, append to existing logs (for ongoing processes)
    """
    # Create logs directory if it doesn't exist
    os.makedirs("logs", exist_ok=True)
    
    if reset:
        # Clear any existing loggers
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
    
    # Mode is "w" if reset is True, otherwise "a"
    file_mode = "w" if reset else "a"
    
    # Configure the detailed log file
    comparison_handler = logging.FileHandler("logs/comparison.log", mode=file_mode)
    comparison_handler.setLevel(logging.DEBUG)
    comparison_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    comparison_handler.setFormatter(comparison_formatter)
    
    # Configure the root logger for higher-level info
    root_handler = logging.FileHandler("logs/comparison-root.log", mode=file_mode)
    root_handler.setLevel(logging.INFO)
    root_formatter = logging.Formatter("%(message)s")
    root_handler.setFormatter(root_formatter)
    
    # Set up console logging
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console_handler.setFormatter(console_formatter)
    
    # Configure the root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    root_logger.addHandler(comparison_handler)
    root_logger.addHandler(root_handler)
    root_logger.addHandler(console_handler)
    
    # Set specific loggers to higher levels to suppress their debug logs
    logging.getLogger("PIL.PngImagePlugin").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
    logging.getLogger("pubchempy").setLevel(logging.WARNING)
    logging.getLogger("watchdog.observers.inotify_buffer").setLevel(logging.WARNING)
    logging.getLogger("watchdog.observers").setLevel(logging.WARNING)
    logging.getLogger("watchdog").setLevel(logging.WARNING)
    
    return root_logger

def get_logger(name):
    """
    Get a logger with the specified name, ensuring the logging system is initialized.
    
    Args:
        name (str): The name of the logger, typically __name__
        
    Returns:
        logging.Logger: A configured logger instance
    """
    # Make sure logs directory exists
    os.makedirs("logs", exist_ok=True)
    
    return logging.getLogger(name)
