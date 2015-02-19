'''
Created on 11/2/2015

@author: victor
'''
import unittest
import numpy
import os
import pyproct.data.handler.featurearray.test.data as data
from pyproct.data.handler.featurearray.featureArrayDataLoader import FeatureArrayDataLoader
from pyproct.data.handler.dataSource import DataSource


class TestFeatArrDataLoader(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.data_path = data.__path__[0]
        
    def compare_loaded_results(self,r1, r2):
        self.assertEqual(r1[1], r2[1])
        d1, d2 = r1[0], r2[0]
        self.assertItemsEqual(d1.keys(), d2.keys())
        for label in d1:
            self.assertItemsEqual(d1[label], d2[label])
    
    def test_load(self):
        loader = FeatureArrayDataLoader(None)
        source = DataSource(os.path.join(TestFeatArrDataLoader.data_path, "numerical_features_1.txt"))
        self.compare_loaded_results(data.expected_1, loader.load_data_from_source(source))
        source.add_info("labels", ["uno", "dos", "tres", "cuatro", "cinco"])
        self.compare_loaded_results(data.expected_2, loader.load_data_from_source(source))
        
    def test_load_one_row(self):
        loader = FeatureArrayDataLoader(None)
        source = DataSource(os.path.join(TestFeatArrDataLoader.data_path, "numerical_features_1.txt"))
        source.add_info("labels", ["uno", "dos", "tres", "cuatro", "cinco"])
        source.add_info("load_only", [2, 3])
        self.compare_loaded_results(data.expected_3, loader.load_data_from_source(source))
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()