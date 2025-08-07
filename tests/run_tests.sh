#!/bin/bash
# Test runner script for Chemical Recommender System

echo "======================================================"
echo "Chemical Recommender System - Test Suite"
echo "======================================================"

# Check if we're in the right directory
if [ ! -f "src/main.py" ]; then
    echo "Error: Please run this script from the project root directory"
    exit 1
fi

# Create logs directory if it doesn't exist
mkdir -p logs

# Set PYTHONPATH to include src directory
export PYTHONPATH="${PYTHONPATH}:$(pwd)/src"

echo "Running Unit Tests..."
echo "------------------------------------------------------"

# Run unit tests with coverage
python -m unittest discover tests/ -p "test_unit_*.py" -v

echo ""
echo "Running Golden Master Tests..."
echo "------------------------------------------------------"
echo "Note: Golden master tests may create baseline files on first run"

# Run golden master tests
python -m unittest tests.test_golden_master -v

