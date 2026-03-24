# Golden master tests - compare against known good outputs
import unittest
import os
import sys
import pandas as pd

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from Comparison.Controller import BatchRun

class TestGoldenMaster(unittest.TestCase):
    """
    Golden Master tests - compare outputs against known good results.
    
    This test validates generated batch output against the committed
    tests/golden.csv baseline.
    """

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
            
        except Exception as e:
            if "Milvus" in str(e) or "connection" in str(e).lower():
                self.skipTest(f"Milvus not available: {e}")
            else:
                raise


if __name__ == '__main__':
    # Instructions for running CSV regression check
    print("\n" + "="*60)
    print("GOLDEN MASTER TESTS")
    print("="*60)
    print("This test compares generated output against tests/golden.csv.")
    print("The golden CSV baseline is read-only during test runs.")
    print("="*60 + "\n")
    
    unittest.main()
