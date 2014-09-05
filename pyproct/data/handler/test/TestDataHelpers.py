"""
Created on Sep 3, 2014

@author: victor
"""
import unittest
from pyproct.data.handler.dataSource import DataSource
from pyproct.data.handler.elementRange import ElementRange
from pyproct.data.handler.sourceGenerator import SourceGenerator
import pyproct.data.handler.test.data as test_data
import os
import glob

class TestDataHelpers(unittest.TestCase):

    def test_data_source_creation(self):
        # Test creation
        source_1 = DataSource("lol.pdb")
        self.assertDictEqual({"source":"lol.pdb"}, source_1.source)
        
        source_2 = DataSource({"source":"lol.pdb", "another_kw": 3})
        self.assertDictEqual({"source":"lol.pdb", "another_kw": 3}, source_2.source)
    
    def test_data_source_comparison(self):
        m_source = DataSource("mini")
        n_source = DataSource("normal")
        another_n_source = DataSource("normal")
        b_source = DataSource("zobig")
        
        # Test lexicographical comparison
        self.assertTrue(m_source < n_source)
        self.assertTrue(n_source == another_n_source)
        self.assertTrue(b_source > n_source)

    def test_element_range_inflation_and_in(self):
        erange = ElementRange(4,9)
        self.assertEqual([4, 5, 6, 7, 8, 9], list(erange))
        
        self.assertTrue(5 in erange)
        self.assertFalse(10 in erange)
        
        c_list = []
        for i in erange:
            c_list.append(i)
        self.assertEqual([4, 5, 6, 7, 8, 9], c_list)
    
    def test_element_range_len(self):
        erange = ElementRange(4,9)
        self.assertEqual(6, len(erange))
        
    def test_source_generator(self):
        """
        TODO: test inflating independently (string, file list, dict) and recursive inflating
        for all.
        """
        class SourceGeneratorMock(SourceGenerator):
            def __init__(self, source_list):
                super(SourceGeneratorMock, self).__init__(source_list)
            
            @classmethod 
            def do_glob(cls, source):
                """
                "Normal" do_glob fails if file or description does not exist. 
                """
                if len(glob.glob(source)) == 0:
                    return [source]
                else:
                    return glob.glob(source)
                    
        sources_list = [
           "file.lel",
           os.path.join(test_data.__path__[0], "lots_of_files.lst"),
           {"source": "protein.pdb"},
           {"source": "another_protein.pdb", "with_kw":"KW!"},
           os.path.join(test_data.__path__[0], "bunch_of_files.lst"),
           {"source": os.path.join(test_data.__path__[0], "*.pdb"), "selection":"all"}
        ]

        expected_ordered_sources = ['file.lel', 
                                    'one.txt', 
                                    {'source': 'two.txt', 'lol': 'lol'}, 
                                    'three.ff', 
                                    {'source': 'protein.pdb'}, 
                                    {'source': 'another_protein.pdb', 'with_kw': 'KW!'}, 
                                    'four.txt', 
                                    'six.txt', 
                                    'seven.txt', 
                                    {'source': 'pdb3.pdb', 'selection': 'all'}, 
                                    {'source': 'pdb2.pdb', 'selection': 'all'}, 
                                    {'source': 'pdb1.pdb', 'selection': 'all'}]

        
        sources =  SourceGeneratorMock.inflate_source_list(sources_list)
        def normalize_sources(sources):
            norm_sources = []
            for s in norm_sources:
                try:
                    s = s['source']
                except Exception:
                    s = s
                norm_sources.append(os.path.basename(s))
            return norm_sources
            
        self.assertItemsEqual(normalize_sources(expected_ordered_sources), 
                              normalize_sources(sources))
                
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()