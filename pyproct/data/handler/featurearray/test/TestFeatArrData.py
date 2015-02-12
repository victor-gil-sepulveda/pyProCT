"""
Created on 11/2/2015

@author: victor
"""
import unittest
import pyproct.data.handler.featurearray.test.data as data
from pyproct.data.handler.featurearray.featureArrayData import FeatureArrayData
import numpy

class TestFeatArrData(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.data_path = data.__path__[0]


    def test_get_all_elements(self):
        expected_2_columns = [[17.0, 91.0], [82.0, 15.0], [26.0, 36.0], [99.0, 51.0], [72.0, 32.0], [35.0, 11.0], [54.0, 38.0], [22.0, 56.0], [20.0, 21.0], [25.0, 34.0], [29.0, 75.0], [94.0, 77.0], [66.0, 98.0], [84.0, 71.0], [55.0, 95.0], [12.0, 4.0], [43.0, 83.0], [1.0, 70.0], [16.0, 33.0]]
        expected_single_column = [[91.0], [15.0], [36.0], [51.0], [32.0], [11.0], [38.0], [56.0], [21.0], [34.0], [75.0], [77.0], [98.0], [71.0], [95.0], [4.0], [83.0], [70.0], [33.0]]

        feature_data = FeatureArrayData(data.expected_3[0])
        numpy.testing.assert_array_equal(expected_2_columns, feature_data.get_all_elements())
        feature_data = FeatureArrayData(data.expected_single_row)
        numpy.testing.assert_array_equal(expected_single_column, feature_data.get_all_elements())
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()