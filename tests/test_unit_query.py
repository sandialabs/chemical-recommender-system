# Unit tests for the query parsing functionality
import unittest
from unittest.mock import patch, Mock, MagicMock
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from Comparison.utils.gen import parseQuery, filterCIDs, fillDict


class TestQueryParsing(unittest.TestCase):
    """Test the query parsing functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.mock_logger = Mock()
    
    @patch('Comparison.utils.gen.pcp.get_properties')
    @patch('Comparison.utils.gen.pcp.get_cids')
    def test_parse_query_by_cid(self, mock_get_cids, mock_get_properties):
        """Test parsing query when input is a CID number"""
        # Mock PubChem responses
        mock_get_properties.return_value = [{
            'SMILES': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'IUPACName': 'aspirin'
        }]
        
        # Test with integer CID
        result = parseQuery(2244)
        
        # Should return [fingerprint, cid, name, smiles]
        self.assertIsNotNone(result[0])  # fingerprint should exist
        self.assertEqual(result[1], 2244)  # CID should match
        self.assertEqual(result[2], 'aspirin')  # name should match
        self.assertEqual(result[3], 'CC(=O)OC1=CC=CC=C1C(=O)O')  # SMILES should match
        
        # Test with string CID
        result = parseQuery("2244")
        self.assertEqual(result[1], 2244)
    
    @patch('Comparison.utils.gen.pcp.get_properties')
    @patch('Comparison.utils.gen.pcp.get_cids')
    def test_parse_query_by_name(self, mock_get_cids, mock_get_properties):
        """Test parsing query when input is a chemical name"""
        # Mock PubChem responses
        mock_get_cids.return_value = [2244]
        mock_get_properties.return_value = [{
            'SMILES': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'IUPACName': 'aspirin'
        }]
        
        result = parseQuery("aspirin")
        
        self.assertIsNotNone(result[0])  # fingerprint should exist
        self.assertEqual(result[1], 2244)
        self.assertEqual(result[2], 'aspirin')
        self.assertEqual(result[3], 'CC(=O)OC1=CC=CC=C1C(=O)O')
    
    @patch('Comparison.utils.gen.pcp.get_compounds')
    def test_parse_query_by_smiles(self, mock_get_compounds):
        """Test parsing query when input is a SMILES string"""
        # Mock PubChem compound object
        mock_compound = Mock()
        mock_compound.cid = 2244
        mock_get_compounds.return_value = [mock_compound]
        
        with patch('Comparison.utils.gen.get_compound_properties') as mock_props:
            mock_props.return_value = ('CC(=O)OC1=CC=CC=C1C(=O)O', 'aspirin')
            
            result = parseQuery("CC(=O)OC1=CC=CC=C1C(=O)O")
            
            self.assertIsNotNone(result[0])
            self.assertEqual(result[1], 2244)
            self.assertEqual(result[2], 'aspirin')
            self.assertEqual(result[3], 'CC(=O)OC1=CC=CC=C1C(=O)O')
    
    def test_parse_query_invalid_input(self):
        """Test parsing with invalid input"""
        with patch('Comparison.utils.gen.pcp.get_cids') as mock_get_cids:
            mock_get_cids.side_effect = Exception("Not found")
            
            with patch('Comparison.utils.gen.pcp.get_compounds') as mock_get_compounds:
                mock_get_compounds.side_effect = Exception("Not found")
                
                result = parseQuery("invalid_input_xyz")
                
                # Should return -1 values for failed parse
                self.assertEqual(result[0], -1)
                self.assertEqual(result[1], -1)
                self.assertEqual(result[2], -1)
                self.assertEqual(result[3], -1)


class TestCIDFiltering(unittest.TestCase):
    """Test the CID filtering functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.sample_found = [
            (0.95, 2244),  # aspirin
            (0.87, 702),   # ethanol  
            (0.72, 297),   # methane
        ]
        
    @patch('Comparison.utils.gen.pcp.get_properties')
    def test_filter_cids_basic(self, mock_get_properties):
        """Test basic CID filtering without restrictions"""
        # Mock PubChem properties response
        mock_get_properties.return_value = [
            {'SMILES': 'CC(=O)OC1=CC=CC=C1C(=O)O', 'IUPACName': 'aspirin', 'MolecularFormula': 'C9H8O4'},
            {'SMILES': 'CCO', 'IUPACName': 'ethanol', 'MolecularFormula': 'C2H6O'},
            {'SMILES': 'C', 'IUPACName': 'methane', 'MolecularFormula': 'CH4'}
        ]
        
        result = filterCIDs(
            found=self.sample_found,
            queryname="test_query",
            incEle=True,  # Include all elements
            smarts_mol=None,  # No substructure filter
            smarts_num=None,
            querysmiles="O",  # Water as query
            job_id="test"
        )
        
        # Should return SMILES for valid compounds
        self.assertEqual(len(result), 3)
        self.assertIn('CC(=O)OC1=CC=CC=C1C(=O)O', result)
        self.assertIn('CCO', result)
        self.assertIn('C', result)
    
    @patch('Comparison.utils.gen.pcp.get_properties')
    def test_filter_cids_element_restriction(self, mock_get_properties):
        """Test CID filtering with element restrictions"""
        mock_get_properties.return_value = [
            {'SMILES': 'CC(=O)OC1=CC=CC=C1C(=O)O', 'IUPACName': 'aspirin', 'MolecularFormula': 'C9H8O4'},
            {'SMILES': 'CCO', 'IUPACName': 'ethanol', 'MolecularFormula': 'C2H6O'},
            {'SMILES': '[Na+].[Cl-]', 'IUPACName': 'sodium chloride', 'MolecularFormula': 'ClNa'}  # Contains Na (not allowed)
        ]
        
        # Add Silicon to allowed elements (Na still not in allowed list)
        result = filterCIDs(
            found=[(0.9, 2244), (0.8, 702), (0.7, 123)],  # Made up CID for sodium chloride
            queryname="test_query", 
            incEle=["Si"],  # Add Silicon to allowed elements (Na not in default or added)
            smarts_mol=None,
            smarts_num=None,
            querysmiles="O",
            job_id="test"
        )
        
        # Sodium chloride should be filtered out (None) because Na is not in allowed elements
        # Default allowed: ["H", "C", "N", "O", "F", "P", "S", "Cl", "Se", "Br", "I"] + ["Si"]
        self.assertEqual(len(result), 3)
        self.assertIn('CC(=O)OC1=CC=CC=C1C(=O)O', result)
        self.assertIn('CCO', result) 
        self.assertIsNone(result[2])  # Sodium chloride should be None (contains Na)


class TestFillDict(unittest.TestCase):
    """Test the dictionary filling functionality"""
    
    @patch('Comparison.utils.gen.filterCIDs')
    def test_fill_dict_basic(self, mock_filter):
        """Test basic dictionary filling"""
        # Mock filter results
        mock_filter.return_value = ['CCO', 'C', None, 'CC(=O)OC1=CC=CC=C1C(=O)O']
        
        comp_dict = {}
        heap = [(0.9, 702), (0.8, 297), (0.7, 999), (0.6, 2244)]  # Last one should be filtered out
        
        result = fillDict(
            comp_dict=comp_dict,
            heap=heap,
            finsize=3,
            queryname="test",
            incEle=True,
            smarts_mol=None,
            smarts_num=None,
            querysmiles="O",
            job_id="test"
        )
        
        # Should have added valid compounds to comp_dict
        self.assertIn(702, comp_dict)  # ethanol
        self.assertIn(297, comp_dict)  # methane
        self.assertIn(2244, comp_dict)  # aspirin
        self.assertNotIn(999, comp_dict)  # filtered out
        
        # Check structure of comp_dict values
        for key, val in comp_dict.items():
            self.assertEqual(len(val), 6)  # [score, score, None, None, None, None]
            self.assertIsInstance(val[0], float)  # similarity score
            self.assertIsInstance(val[1], float)  # structural score
        
        # Should return smiles_dict
        self.assertIsInstance(result, dict)
        self.assertIn(702, result)
        self.assertEqual(result[702], 'CCO')


if __name__ == '__main__':
    unittest.main()
