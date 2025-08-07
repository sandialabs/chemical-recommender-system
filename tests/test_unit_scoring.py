# Unit tests for comparison helpers (scoring functions)
import unittest
from unittest.mock import patch, Mock, MagicMock
import pandas as pd
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from Comparison.Helpers.StrucComp import extraStrucComp
from Comparison.Helpers.ThermComp import thermalComparison
from Comparison.Helpers.ToxicComp import toxicComparison


class TestStructuralComparison(unittest.TestCase):
    """Test structural comparison scoring"""
    
    def setUp(self):
        """Set up test data"""
        self.comp_dict = {
            '2244': [0.9, 0.9, None, None, None, None],  # aspirin (query)
            '702': [0.8, 0.8, None, None, None, None],   # ethanol
            '297': [0.7, 0.7, None, None, None, None],   # methane
        }
        
        # Mock CSV data
        self.mock_df = pd.DataFrame({
            'MoleculeID': [2244, 702, 297],
            'MolWeight': [180.16, 46.07, 16.04],
            'nbRing': [2, 0, 0],
            'nbLipinskiFailures': [0, 0, 0], 
            'TopoPolSurfAir': [63.6, 20.2, 0.0],
            'nbC': [9, 2, 1]
        })
    
    @patch('Comparison.Helpers.StrucComp.pd.read_csv')
    def test_structural_comparison_basic(self, mock_read_csv):
        """Test basic structural comparison calculation"""
        mock_read_csv.return_value = self.mock_df
        
        result = extraStrucComp(self.comp_dict.copy(), job_id="test")
        
        # Query should have molecular weight score of 1
        self.assertEqual(result['2244'][2], 1)
        
        # Other compounds should have calculated molecular weight scores
        self.assertIsNotNone(result['702'][2])
        self.assertIsNotNone(result['297'][2])
        
        # Scores should be between 0 and 1
        for cid, values in result.items():
            if values[2] is not None:
                self.assertGreaterEqual(values[2], 0)
                self.assertLessEqual(values[2], 1)
    
    @patch('Comparison.Helpers.StrucComp.pd.read_csv')
    def test_structural_comparison_missing_columns(self, mock_read_csv):
        """Test structural comparison with missing columns"""
        # Create DataFrame without some structural columns
        incomplete_df = pd.DataFrame({
            'MoleculeID': [2244, 702, 297],
            'MolWeight': [180.16, 46.07, 16.04],
            # Missing: nbRing, nbLipinskiFailures, TopoPolSurfAir, nbC
        })
        mock_read_csv.return_value = incomplete_df
        
        # Should not raise error, just skip missing comparisons
        result = extraStrucComp(self.comp_dict.copy(), job_id="test")
        
        # Should still have molecular weight scores
        self.assertEqual(result['2244'][2], 1)
        self.assertIsNotNone(result['702'][2])


class TestThermalComparison(unittest.TestCase):
    """Test thermal property comparison scoring"""
    
    def setUp(self):
        """Set up test data"""
        self.comp_dict = {
            '2244': [0.9, 0.9, 0.9, None, None, None],  # aspirin
            '702': [0.8, 0.8, 0.8, None, None, None],   # ethanol
            '297': [0.7, 0.7, 0.7, None, None, None],   # methane
        }
        
        self.mock_df = pd.DataFrame({
            'MoleculeID': [1234, 2244, 702, 297],  # Query CID first, then candidates
            'MP_pred': [100.0, 135.0, -114.0, -182.0],  # Query value first
            'BP_pred': [150.0, 140.0, 78.0, -161.0], 
            'LogP_pred': [1.0, 1.19, -0.31, 1.09],
            'LogHL_pred': [-3.0, -5.2, -1.8, 1.2],
            'LogVP_pred': [0.0, -2.3, 1.2, 3.8]
        })
    
    @patch('Comparison.Helpers.ThermComp.pd.read_csv')
    def test_thermal_comparison_all_properties(self, mock_read_csv):
        """Test thermal comparison with all properties enabled"""
        mock_read_csv.return_value = self.mock_df
        
        # Enable all thermal properties
        tarray = [True, True, True, True, True]  # MP, BP, LogP, Hlaw, VP
        
        result = thermalComparison(self.comp_dict.copy(), tarray, job_id="test")
        
        # All compounds should have thermal scores
        for cid, values in result.items():
            self.assertIsNotNone(values[3])  # thermal score at index 3
            self.assertGreaterEqual(values[3], 0)
            self.assertLessEqual(values[3], 1)
    
    @patch('Comparison.Helpers.ThermComp.pd.read_csv')
    def test_thermal_comparison_selective_properties(self, mock_read_csv):
        """Test thermal comparison with only some properties enabled"""
        mock_read_csv.return_value = self.mock_df
        
        # Enable only MP and LogP
        tarray = [True, False, True, False, False]
        
        result = thermalComparison(self.comp_dict.copy(), tarray, job_id="test")
        
        # Should still calculate thermal scores
        for cid, values in result.items():
            self.assertIsNotNone(values[3])
    
    @patch('Comparison.Helpers.ThermComp.pd.read_csv')
    def test_thermal_comparison_no_properties(self, mock_read_csv):
        """Test thermal comparison with no properties enabled"""
        mock_read_csv.return_value = self.mock_df
        
        # Disable all thermal properties
        tarray = [False, False, False, False, False]
        
        result = thermalComparison(self.comp_dict.copy(), tarray, job_id="test")
        
        # Should still have thermal scores (all 1.0 since no comparisons)
        for cid, values in result.items():
            self.assertEqual(values[3], 1.0)


class TestToxicityComparison(unittest.TestCase):
    """Test toxicity comparison scoring"""
    
    def setUp(self):
        """Set up test data"""
        self.comp_dict = {
            '2244': [0.9, 0.9, 0.9, 0.9, None, None],  # aspirin
            '702': [0.8, 0.8, 0.8, 0.8, None, None],   # ethanol
            '297': [0.7, 0.7, 0.7, 0.7, None, None],   # methane
        }
        
        self.mock_df = pd.DataFrame({
            'MoleculeID': [2244, 702, 297],
            'LogBCF_pred': [0.8, -0.5, 0.3],
            'CATMoS_EPA_pred': [0.2, 0.1, 0.05],
            'CATMoS_LD50_pred': [500.0, 1000.0, 800.0]
        })
    
    @patch('Comparison.Helpers.ToxicComp.pd.read_csv')
    def test_toxicity_comparison_basic(self, mock_read_csv):
        """Test basic toxicity comparison calculation"""
        mock_read_csv.return_value = self.mock_df
        
        result = toxicComparison(self.comp_dict.copy(), job_id="test")
        
        # All compounds should have toxicity scores
        for cid, values in result.items():
            self.assertIsNotNone(values[4])  # toxicity score at index 4
            self.assertIsInstance(values[4], (int, float))  # Should be numeric
    
    @patch('Comparison.Helpers.ToxicComp.pd.read_csv')
    def test_toxicity_comparison_missing_data(self, mock_read_csv):
        """Test toxicity comparison with missing data"""
        # DataFrame with some missing values
        incomplete_df = pd.DataFrame({
            'MoleculeID': [2244, 702, 297],
            'LogBCF_pred': [0.8, "N/A", 0.3],  # Missing BCF for ethanol
            'CATMoS_EPA_pred': [0.2, 0.1, 0.05],
            'CATMoS_LD50_pred': [500.0, 1000.0, 800.0]
        })
        mock_read_csv.return_value = incomplete_df
        
        result = toxicComparison(self.comp_dict.copy(), job_id="test")
        
        # Should use default values for missing data
        self.assertIsNotNone(result['702'][4])  # Should have a score despite missing BCF
    
    @patch('Comparison.Helpers.ToxicComp.pd.read_csv')
    def test_toxicity_comparison_missing_columns(self, mock_read_csv):
        """Test toxicity comparison with missing columns"""
        # DataFrame missing some toxicity columns
        incomplete_df = pd.DataFrame({
            'MoleculeID': [2244, 702, 297],
            'LogBCF_pred': [0.8, -0.5, 0.3],
            # Missing: CATMoS_EPA_pred, CATMoS_LD50_pred
        })
        mock_read_csv.return_value = incomplete_df
        
        result = toxicComparison(self.comp_dict.copy(), job_id="test")
        
        # Should use default values for missing columns
        for cid, values in result.items():
            self.assertIsNotNone(values[4])


class TestIntegratedScoring(unittest.TestCase):
    """Test the integration of all scoring components"""
    
    def setUp(self):
        """Set up comprehensive test data"""
        self.comp_dict = {
            '1234': [1.0, 1.0, None, None, None, None],  # query compound
            '2244': [0.9, 0.9, None, None, None, None],  # aspirin
            '702': [0.8, 0.8, None, None, None, None],   # ethanol
        }
        
        self.mock_df = pd.DataFrame({
            'MoleculeID': [1234, 2244, 702],  # Query first, then candidates
            'MolWeight': [150.0, 180.16, 46.07],
            'MP_pred': [100.0, 135.0, -114.0],  # Query value first
            'BP_pred': [140.0, 140.0, 78.0],
            'LogP_pred': [1.0, 1.19, -0.31],
            'LogBCF_pred': [0.5, 0.8, -0.5],
            'CATMoS_EPA_pred': [0.15, 0.2, 0.1],
            'CATMoS_LD50_pred': [750.0, 500.0, 1000.0],
            'nbRing': [1, 2, 0],
            'nbLipinskiFailures': [0, 0, 0],
            'TopoPolSurfAir': [40.0, 63.6, 20.2],
            'nbC': [6, 9, 2]
        })
    
    @patch('Comparison.Helpers.ToxicComp.pd.read_csv')
    @patch('Comparison.Helpers.ThermComp.pd.read_csv')
    @patch('Comparison.Helpers.StrucComp.pd.read_csv')
    def test_complete_scoring_pipeline(self, mock_struc_csv, mock_therm_csv, mock_toxic_csv):
        """Test complete scoring pipeline"""
        # Mock all CSV reads
        mock_struc_csv.return_value = self.mock_df
        mock_therm_csv.return_value = self.mock_df
        mock_toxic_csv.return_value = self.mock_df
        
        comp_dict = self.comp_dict.copy()
        
        # Run all scoring functions
        comp_dict = extraStrucComp(comp_dict, job_id="test")
        comp_dict = thermalComparison(comp_dict, [True, True, True, True, True], job_id="test")
        comp_dict = toxicComparison(comp_dict, job_id="test")
        
        # Verify all scores are calculated
        for cid, values in comp_dict.items():
            self.assertEqual(len(values), 6)
            
            if cid == '1234':  # Query CID - only check molecular weight
                self.assertIsNotNone(values[2])  # molecular weight should be 1
                # Query doesn't get thermal or toxicity scores compared to itself
            else:  # Candidate CIDs - check all scores
                self.assertIsNotNone(values[2])  # molecular weight
                self.assertIsNotNone(values[3])  # thermal
                self.assertIsNotNone(values[4])  # toxicity
            # values[5] will be filled by SA score later


if __name__ == '__main__':
    unittest.main()
