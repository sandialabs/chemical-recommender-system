# Golden master tests - compare against known good outputs
import unittest
import json
import os
import sys
from unittest.mock import patch, Mock
import pandas as pd

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from Comparison.Controller import BatchRun

class TestGoldenMaster(unittest.TestCase):
    """
    Golden Master tests - compare outputs against known good results.
    
    These tests use pre-recorded 'golden' outputs to ensure the system
    produces consistent results. When you make changes that intentionally
    change outputs, you'll need to update the golden files.
    """
    
    def setUp(self):
        """Set up paths for golden master files"""
        self.test_data_dir = os.path.join(os.path.dirname(__file__), 'golden_data')
        os.makedirs(self.test_data_dir, exist_ok=True)
    
    def save_golden_master(self, test_name, data):
        """Save golden master data to file"""
        filepath = os.path.join(self.test_data_dir, f"{test_name}.json")
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2, default=str)
    
    def load_golden_master(self, test_name):
        """Load golden master data from file"""
        filepath = os.path.join(self.test_data_dir, f"{test_name}.json")
        if not os.path.exists(filepath):
            return None
        with open(filepath, 'r') as f:
            return json.load(f)
    
    def assert_golden_master(self, test_name, actual_data, tolerance=1e-6):
        """Compare actual data against golden master"""
        expected_data = self.load_golden_master(test_name)
        
        if expected_data is None:
            # First time running - save as golden master
            self.save_golden_master(test_name, actual_data)
            print(f"Saved new golden master: {test_name}")
            return
        
        # Compare the data structures
        self.compare_data_structures(expected_data, actual_data, tolerance)
    
    def compare_data_structures(self, expected, actual, tolerance):
        """Recursively compare data structures with tolerance for floats"""
        if isinstance(expected, dict) and isinstance(actual, dict):
            self.assertEqual(set(expected.keys()), set(actual.keys()))
            for key in expected.keys():
                self.compare_data_structures(expected[key], actual[key], tolerance)
        elif isinstance(expected, list) and isinstance(actual, list):
            self.assertEqual(len(expected), len(actual))
            for exp_item, act_item in zip(expected, actual):
                self.compare_data_structures(exp_item, act_item, tolerance)
        elif isinstance(expected, float) and isinstance(actual, float):
            self.assertAlmostEqual(expected, actual, delta=tolerance)
        else:
            self.assertEqual(expected, actual)

    def test_exact_cid_6517_batch_reproduction(self):
        """
        Test that reproduces your exact CID 6517 search using BatchRun
        Parameters: 6517, 30, [True,True,False,False,False], False, [Si], 1, CCO, 1
        """
        # Skip if no golden CSV available
        golden_file = os.path.join(os.path.dirname(__file__), 'golden.csv')
        if not os.path.exists(golden_file):
            self.skipTest("golden.csv not found - cannot validate against known results")
        
        try:
            # Create the exact batch input from your parameters
            # Format: query, final_number, thermo_array, include_all_elements, include_specific_elements, disallow_isotopes, substructure_search, number_substructure_search
            batch_input = "6517, 30, [True,True,False,False,False], False, [Si], 1, CCO, 1"
            
            # Run BatchRun - it processes the query and saves results to CSV files
            BatchRun(
                batch_text=batch_input,
                containers=[],  # No additional containers
                job_id="test_golden_cid_6517_batch"
            )
            
            # BatchRun generates data.csv from the last SingleRun
            # Read the generated CSV file
            generated_csv = "src/Comparison/LocalIO/data.csv"
            if not os.path.exists(generated_csv):
                self.fail("BatchRun did not generate expected data.csv file")
            
            generated_df = pd.read_csv(generated_csv)
            golden_df = pd.read_csv(golden_file)
            
            # Validate basic structure
            self.assertEqual(len(generated_df), len(golden_df), 
                           f"Expected {len(golden_df)} results, got {len(generated_df)}")
            
            # The first result should be the query (CID 6517)
            query_result = generated_df.iloc[0]
            self.assertEqual(int(query_result['CID']), 6517, "First result should be query CID 6517")
            
            # Query should have None for overall score (no self-comparison)
            self.assertTrue(pd.isna(query_result['Overall']), "Query should have None overall score")
            
            # Validate the ranking matches your golden data exactly
            actual_ranked_results = generated_df[generated_df['CID'] != 6517].copy()  # Skip query
            golden_ranked = golden_df[golden_df['CID'] != 6517].copy()  # Skip query from golden
            
            # Check that we have the same results in the same order
            for i, (actual_row, golden_row) in enumerate(zip(actual_ranked_results.itertuples(), golden_ranked.itertuples())):
                actual_cid = int(actual_row.CID)
                golden_cid = int(golden_row.CID)
                
                self.assertEqual(actual_cid, golden_cid, 
                               f"Result #{i+1}: Expected CID {golden_cid}, got {actual_cid}")
                
                # Check overall score (with tolerance for floating point differences)
                if not pd.isna(actual_row.Overall) and not pd.isna(golden_row.Overall):
                    self.assertAlmostEqual(float(actual_row.Overall), float(golden_row.Overall), 
                                         places=5,  # Slightly less strict for batch processing
                                         msg=f"CID {actual_cid}: Overall score mismatch")
                
                # Check individual component scores
                score_names = ['Fingerprint', 'Molecular', 'Thermophysical', 'Toxicity', 'Synthetic']
                
                for score_name in score_names:
                    if (hasattr(actual_row, score_name) and hasattr(golden_row, score_name) and
                        not pd.isna(getattr(actual_row, score_name)) and not pd.isna(getattr(golden_row, score_name))):
                        actual_score = float(getattr(actual_row, score_name))
                        golden_score = float(getattr(golden_row, score_name))
                        self.assertAlmostEqual(actual_score, golden_score, places=4,  
                                             msg=f"CID {actual_cid}: {score_name} score mismatch")
            
            # Validate specific expected results from your golden.csv
            top_5_cids = list(actual_ranked_results['CID'].head(5).astype(int))
            expected_top_5 = [524710, 20663005, 21515778, 19877451, 59678459]  # From your golden.csv
            
            self.assertEqual(top_5_cids, expected_top_5, 
                           "Top 5 results should match golden data exactly")
            
            # Validate top result details
            top_result = actual_ranked_results.iloc[0]
            golden_top = golden_df[golden_df['CID'] == 524710].iloc[0]
            
            self.assertAlmostEqual(float(top_result['Overall']), float(golden_top['Overall']), places=5,
                                 msg="Top result overall score should match golden data")
            self.assertAlmostEqual(float(top_result['Fingerprint']), float(golden_top['Fingerprint']), places=4,
                                 msg="Top result fingerprint score should match golden data")
            
            # Store comprehensive summary for regression testing
            test_summary = {
                'query_cid': 6517,
                'batch_parameters': {
                    'finnum': 30,
                    'tarray': [True, True, False, False, False],
                    'incEle': False,
                    'include_specific_elements': ['Si'],
                    'disallow_isotopes': 1,
                    'substructure': 'CCO',
                    'substructure_count': 1
                },
                'num_results': len(actual_ranked_results),
                'all_result_cids': list(actual_ranked_results['CID'].astype(int)),
                'all_overall_scores': [float(x) if not pd.isna(x) else None 
                                     for x in actual_ranked_results['Overall']],
                'top_10_details': [
                    {
                        'cid': int(row.CID),
                        'overall': float(row.Overall) if not pd.isna(row.Overall) else None,
                        'fingerprint': float(row.Fingerprint) if hasattr(row, 'Fingerprint') and not pd.isna(row.Fingerprint) else None,
                        'molecular': float(row.Molecular) if hasattr(row, 'Molecular') and not pd.isna(row.Molecular) else None,
                        'thermal': float(row.Thermophysical) if hasattr(row, 'Thermophysical') and not pd.isna(row.Thermophysical) else None,
                        'toxicity': float(row.Toxicity) if hasattr(row, 'Toxicity') and not pd.isna(row.Toxicity) else None,
                        'synthetic': float(row.Synthetic) if hasattr(row, 'Synthetic') and not pd.isna(row.Synthetic) else None
                    }
                    for row in actual_ranked_results.head(10).itertuples()
                ],
                'validation_status': 'exact_match_with_golden_csv',
                'test_timestamp': str(pd.Timestamp.now())
            }
            
            self.assert_golden_master("exact_cid_6517_batch_reproduction", test_summary, tolerance=1e-5)
                
        except Exception as e:
            if "Milvus" in str(e) or "connection" in str(e).lower():
                self.skipTest(f"Milvus not available: {e}")
            else:
                raise


if __name__ == '__main__':
    # Instructions for updating golden masters
    print("\n" + "="*60)
    print("GOLDEN MASTER TESTS")
    print("="*60)
    print("These tests compare outputs against saved 'golden' results.")
    print("If this is your first time running, golden masters will be created.")
    print("If you've made intentional changes, delete the golden_data/ folder")
    print("to regenerate golden masters.")
    print("="*60 + "\n")
    
    unittest.main()
