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

    def test_load_one(self):
        loader = FeatureArrayDataLoader(None)
        source = DataSource(os.path.join(TestFeatArrDataLoader.data_path, "numerical_features_1.txt"))
        self.assertSequenceEqual(data.expected_1, loader.load_data_from_source(source))
        source.add_info("labels", ["uno", "dos", "tres", "cuatro", "cinco"])
        self.assertSequenceEqual(data.expected_2, loader.load_data_from_source(source))
        
    def test_load_one_one_row(self):
        loader = FeatureArrayDataLoader(None)
        source = DataSource(os.path.join(TestFeatArrDataLoader.data_path, "numerical_features_1.txt"))
        source.add_info("labels", ["uno", "dos", "tres", "cuatro", "cinco"])
        source.add_info("load_only", [2, 3])
        self.assertSequenceEqual(data.expected_3, loader.load_data_from_source(source))
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()