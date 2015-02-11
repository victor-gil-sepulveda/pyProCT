"""
Created on 11/2/2015

@author: victor
"""
import unittest
import pyproct.data.handler.featurearray.test.data as data
from pyproct.data.handler.featurearray.featureArrayData import FeatureArrayData
from scipy.spatial.distance import  pdist
class TestFeatArrData(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.data_path = data.__path__[0]


    def test_get_all_elements(self):
        feature_data = FeatureArrayData(data.expected_3[0])
        print feature_data.get_all_elements()
        print pdist(feature_data.get_all_elements())
        feature_data = FeatureArrayData(data.expected_single_row)
        print feature_data.get_all_elements()
        print pdist(feature_data.get_all_elements())
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()