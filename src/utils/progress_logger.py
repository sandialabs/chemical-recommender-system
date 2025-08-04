# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

import logging
import json
import threading
import os
from typing import Optional, Dict, Any


class ProgressLogger:
    """
    A unified logger that handles both traditional file logging and real-time progress updates.
    
    This logger automatically sends appropriate messages to both log files and progress queues
    for live UI updates, eliminating the need for separate logging and progress systems.
    """
    
    # Class-level storage for logger instances and progress queues
    _instances: Dict[str, 'ProgressLogger'] = {}
    _progress_queues: Dict[str, Any] = {}
    _lock = threading.Lock()
    
    def __init__(self, job_id: str, progress_queue=None):
        self.job_id = job_id
        
        # Store progress queue for this job_id so other calls can access it
        if progress_queue:
            ProgressLogger._progress_queues[job_id] = progress_queue
            self.progress_queue = progress_queue
        elif job_id in ProgressLogger._progress_queues:
            # Use existing progress queue for this job_id
            self.progress_queue = ProgressLogger._progress_queues[job_id]
        else:
            self.progress_queue = None
        
        # Create a standard logger for this job
        self.logger = logging.getLogger(f"progress.{job_id}")
        
        # Ensure logs directory exists
        os.makedirs("logs", exist_ok=True)
        
        # Set up file handlers if not already configured
        if not self.logger.handlers:
            # Detailed log file
            detail_handler = logging.FileHandler(f"logs/job_{job_id}.log", mode="w")
            detail_handler.setLevel(logging.DEBUG)
            detail_formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
            detail_handler.setFormatter(detail_formatter)
            
            # Also add to main comparison log
            main_handler = logging.FileHandler("logs/comparison.log", mode="a")
            main_handler.setLevel(logging.INFO)
            main_formatter = logging.Formatter(
                f"[{job_id}] %(asctime)s - %(levelname)s - %(message)s"
            )
            main_handler.setFormatter(main_formatter)
            
            self.logger.addHandler(detail_handler)
            self.logger.addHandler(main_handler)
            self.logger.setLevel(logging.DEBUG)
    
    @classmethod
    def get_logger(cls, job_id: str, progress_queue=None) -> 'ProgressLogger':
        """
        Get or create a ProgressLogger instance for the given job_id.
        
        Args:
            job_id: Unique identifier for the job/search session
            progress_queue: Optional queue for real-time progress updates
            
        Returns:
            ProgressLogger instance for the job
        """
        with cls._lock:
            # Store progress queue for this job_id if provided
            if progress_queue:
                cls._progress_queues[job_id] = progress_queue
            
            if job_id not in cls._instances:
                cls._instances[job_id] = cls(job_id, progress_queue)
            else:
                # Update existing instance with progress queue if available
                existing_instance = cls._instances[job_id]
                if job_id in cls._progress_queues and not existing_instance.progress_queue:
                    existing_instance.progress_queue = cls._progress_queues[job_id]
            
            return cls._instances[job_id]
    
    @classmethod
    def cleanup_job(cls, job_id: str):
        """Clean up resources for a completed job."""
        with cls._lock:
            if job_id in cls._instances:
                # Remove handlers to prevent memory leaks
                logger_instance = cls._instances[job_id]
                for handler in logger_instance.logger.handlers[:]:
                    handler.close()
                    logger_instance.logger.removeHandler(handler)
                
                del cls._instances[job_id]
            
            if job_id in cls._progress_queues:
                del cls._progress_queues[job_id]
    
    def _send_progress(self, status: str, message: str):
        """Send progress update to the UI if queue is available."""
        if self.progress_queue:
            try:
                progress_msg = json.dumps({
                    "status": status,
                    "detail": message
                })
                self.progress_queue.put_nowait(progress_msg)
                self.logger.debug(f"Progress sent for job {self.job_id}: {status} - {message}")
            except Exception as e:
                # Don't break logging if progress fails
                self.logger.debug(f"Failed to send progress update: {e}")
        else:
            self.logger.debug(f"No progress queue available for job {self.job_id} - message: {status} - {message}")
    
    def info(self, message: str, progress_status: Optional[str] = None):
        """
        Log an info message and optionally send progress update.
        
        Args:
            message: The log message
            progress_status: Optional status for UI progress (e.g., "Running", "Processing")
        """
        self.logger.info(message)
        if progress_status:
            self._send_progress(progress_status, message)
    
    def debug(self, message: str):
        """Log a debug message (debug messages don't go to progress)."""
        self.logger.debug(message)
    
    def warning(self, message: str, progress_status: Optional[str] = None):
        """
        Log a warning and optionally send progress update.
        
        Args:
            message: The warning message
            progress_status: Optional status for UI progress (typically "Warning")
        """
        self.logger.warning(message)
        if progress_status:
            self._send_progress(progress_status, message)
    
    def error(self, message: str, send_to_progress: bool = True):
        """
        Log an error and optionally send to progress UI.
        
        Args:
            message: The error message
            send_to_progress: Whether to send this error to the UI progress
        """
        self.logger.error(message)
        if send_to_progress:
            self._send_progress("Error", message)
    
    def progress(self, status: str, message: str):
        """
        Send a progress-only update (logs at info level but always sends to UI).
        
        Args:
            status: Progress status (e.g., "Running", "Processing", "Finished")
            message: Progress message for the UI
        """
        self.logger.info(f"[PROGRESS] {status}: {message}")
        self._send_progress(status, message)
    
    def step(self, step_name: str, status: str = "Running"):
        """
        Log the start of a processing step.
        
        Args:
            step_name: Name of the step being started
            status: Progress status (default: "Running")
        """
        message = f"Starting: {step_name}"
        self.info(message, status)
    
    def step_complete(self, step_name: str, status: str = "Processing"):
        """
        Log the completion of a processing step.
        
        Args:
            step_name: Name of the step that completed
            status: Progress status (default: "Processing")
        """
        message = f"Completed: {step_name}"
        self.info(message, status)
    
    def progress_update(self, percentage: float, activity: str, status: str = "Processing"):
        """
        Send a percentage-based progress update.
        
        Args:
            percentage: Completion percentage (0-100)
            activity: What is currently being processed
            status: Progress status (default: "Processing")
        """
        message = f"{activity} ({percentage:.1f}% complete)"
        self.info(message, status)


# Convenience function for backward compatibility
def get_progress_logger(job_id: str, progress_queue=None) -> ProgressLogger:
    """
    Convenience function to get a ProgressLogger instance.
    
    Args:
        job_id: Unique identifier for the job/search session
        progress_queue: Optional queue for real-time progress updates
        
    Returns:
        ProgressLogger instance for the job
    """
    return ProgressLogger.get_logger(job_id, progress_queue)


# Context manager for handling job lifecycle
class ProgressContext:
    """Context manager for automatic job cleanup."""
    
    def __init__(self, job_id: str, progress_queue=None):
        self.job_id = job_id
        self.logger = ProgressLogger.get_logger(job_id, progress_queue)
    
    def __enter__(self) -> ProgressLogger:
        return self.logger
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type:
            self.logger.error(f"Job failed with error: {exc_val}")
        ProgressLogger.cleanup_job(self.job_id)


def progress_context(job_id: str, progress_queue=None):
    """
    Create a context manager for progress logging with automatic cleanup.
    
    Usage:
        with progress_context("job_123", queue) as logger:
            logger.step("Processing data")
            # ... do work ...
            logger.step_complete("Processing data")
    """
    return ProgressContext(job_id, progress_queue)
