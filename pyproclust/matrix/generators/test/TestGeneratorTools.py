'''
Created on 16/03/2012

@author: victor
'''
import unittest
from pyproclust.matrix.generators.generatorTools import get_nosymm_ranges,\
                                                get_symm_ranges

class Test(unittest.TestCase):

    def test_get_nosymm_ranges(self):
        expected_ranges = [(0, 2), (2, 2), (4, 3), (7, 4), (11, 4)] 
        self.assertEqual(expected_ranges,get_nosymm_ranges(5,16))
        
        expected_ranges = [(0, 4), (4, 5), (9, 9)]
        self.assertEqual(expected_ranges,get_nosymm_ranges(3,19))
        
    def test_get_symm_ranges(self):
        expected_ranges = [(0, 3), (3, 3), (6, 3), (9, 3), (12, 4)]
        self.assertEqual(expected_ranges,get_symm_ranges(5,16))
        
        expected_ranges = [(0, 2), (2, 2), (4, 2), (6, 2), (8, 2), (10,2), (12,2)]
        self.assertEqual(expected_ranges,get_symm_ranges(8,14))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()