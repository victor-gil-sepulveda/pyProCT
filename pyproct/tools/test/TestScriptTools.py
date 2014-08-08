"""
Created on 02/05/2012

@author: victor
"""
import unittest
from pyproct.tools.scriptTools import get_not_repeated_file_name,\
    create_directory
import os
import  pyproct.tools.test.data as data_mod

class testScriptTools(unittest.TestCase):
    
    def test_get_not_repeated_file_name(self):
        path =  data_mod.__path__[0]
        existing_file = os.path.join(path, "pdb_data.py")
        re_existing_file = os.path.join(path, "_pdb_data.py")
        non_existing_file = os.path.join(path, "nonexisting.txt")
        self.assertEquals(non_existing_file, get_not_repeated_file_name(non_existing_file))
        self.assertEquals(re_existing_file, get_not_repeated_file_name(existing_file))
    
    def test_create_directory(self):
        create_directory("tmp_test/test")
        self.assertTrue(os.path.exists("tmp_test/test"))
        os.system("rm -rf tmp_test")
        self.assertFalse(create_directory("/folder_at_root", True))
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_separate_pure_from_mixed_clusters']
    unittest.main()