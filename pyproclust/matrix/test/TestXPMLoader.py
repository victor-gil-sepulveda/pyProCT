'''
Created on 15/03/2012

@author: victor
'''
import unittest

from pyproclust.matrix import xpmLoader
import pyproclust.matrix.generators.test.data as test_data
import cStringIO
from pyproclust.matrix.completeMatrix import CompleteDistanceMatrix

class Test(unittest.TestCase):

    def test_load_and_convert(self):
        input_handler = cStringIO.StringIO(test_data.xpm_matrix_string)
        loader = xpmLoader.XPMConverter()
        complete_matrix = loader.convert(input_handler)
        expected_matrix = CompleteDistanceMatrix(test_data.xpm_values)
        self.assertEqual(True, complete_matrix == expected_matrix)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_load_and_convert']
    unittest.main()