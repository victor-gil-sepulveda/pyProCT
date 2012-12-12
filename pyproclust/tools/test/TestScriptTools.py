'''
Created on 02/05/2012

@author: victor
'''
import unittest
from pyproclust.tools.scriptTools import get_not_repeated_file_name


class testScriptTools(unittest.TestCase):
    
    def test_get_not_repeated_file_name(self):
        existing_file = "data/pdb_data.py"
        re_existing_file = "data/_pdb_data.py"
        non_existing_file = "data/nonexisting.txt"
        self.assertEquals(non_existing_file, get_not_repeated_file_name(non_existing_file))
        self.assertEquals(re_existing_file, get_not_repeated_file_name(existing_file))
    
    def test_create_directory(self):
        self.fail("TODO")
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_separate_pure_from_mixed_clusters']
    unittest.main()