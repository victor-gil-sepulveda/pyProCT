'''
Created on 5/9/2014

@author: victor
'''
import unittest
import pyproct.data.handler.protein.test.data as test_data
import os
from pyproct.data.handler.protein.proteinEnsembleDataLoader import ProteinEnsembleDataLoader
from pyproct.data.handler.dataSource import DataSource


class TestProteinEnsembleDataLoader(unittest.TestCase):


    def test_load(self):
        loader = ProteinEnsembleDataLoader()
        source1 = DataSource(os.path.join(test_data.__path__[0], "pdb1.pdb"))
        print list(loader.load(source1))
        source2 = DataSource(os.path.join(test_data.__path__[0], "pdb2.pdb"))
        print list(loader.load(source2))
        source3 = DataSource(os.path.join(test_data.__path__[0], "pdb3.pdb"))
        print list(loader.load(source3))
        
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()