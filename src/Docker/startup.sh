#!/bin/bash

# Set application root directory for consistent paths
export APP_ROOT=/app
export FLASK_DEBUG=1
export FLASK_APP=src/App/App.py
export PYTHONPATH="${APP_ROOT}:${PYTHONPATH}"

# Start the Flask app using Gunicorn with sync worker (gevent causes issues with Milvus client)
# Timeout set to 0 to disable worker timeouts (needed for long-running chemical searches)
exec gunicorn src.App.App:app --worker-class sync --bind 0.0.0.0:5005 --timeout 0 --workers 1
