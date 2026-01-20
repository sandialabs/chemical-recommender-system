# Test README - How to run and understand the tests

## Chemical Recommender System Tests

This directory contains tests for the Chemical Recommender System (CRS). The tests are designed to be simple but comprehensive.

### Test Types

#### 1. Unit Tests (`test_unit_*.py`)
Test individual functions in isolation using mocks:
- **`test_unit_query.py`**: Tests query parsing, CID filtering, and dictionary filling
- **`test_unit_scoring.py`**: Tests structural, thermal, and toxicity scoring functions

#### 2. Golden Master Tests (`test_golden_master.py`)
Compare current outputs against saved "golden" reference outputs:
- Ensures consistent behavior over time
- Catches unintended changes to results
- Creates baseline files on first run

### Running Tests

#### Quick Start
```bash
# From project root directory
cd /path/to/chemical-recommender-system
./tests/run_tests.sh
```

#### Run without Docker (local Python, no Golden Test)
- Clone the repo and enter it: `git clone https://github.com/sandialabs/chemical-recommender-system.git CRS && cd CRS`
- Create a venv and install minimal deps: `python3 -m venv .venv && source .venv/bin/activate && pip install pytest pandas pubchempy rdkit`
- Run focused tests (no Milvus/OPERA needed): `pytest tests/test_unit_query.py` and `pytest tests/test_unit_scoring.py`

#### Individual Test Files
```bash
# Run specific test files
python -m unittest tests.test_unit_query -v
python -m unittest tests.test_unit_scoring -v
python -m unittest tests.test_golden_master -v
```

#### Unit Tests Verify:
- ✅ Query parsing (CID, name, SMILES)
- ✅ Candidate filtering (elements, substructures)
- ✅ Scoring calculations (structural, thermal, toxicity)
- ✅ Error handling (missing data, invalid inputs)

#### Golden Master Tests Verify:
- ✅ Consistent ranking of compounds
- ✅ Stable numerical outputs
- ✅ Complete pipeline behavior
- ✅ Integration between components

### Adding New Tests

#### For a New Function:
1. Add unit test in appropriate `test_unit_*.py` file
2. Mock external dependencies
3. Test normal cases and edge cases

#### For End-to-End Behavior:
1. Add golden master test in `test_golden_master.py`
2. Run once to establish baseline
3. Future runs will verify consistency
