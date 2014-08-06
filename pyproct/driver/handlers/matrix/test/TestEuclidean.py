"""
Created on 28/03/2013

@author: victor
"""
import unittest
from pyproct.driver.handlers.matrix.euclideanDistanceMatrixBuilder import EuclideanDistanceMatrixBuilder
import numpy

class TestEuclidean(unittest.TestCase):


    def test_calculate_geom_center(self):
        coordsets = numpy.array([
                     [
                      [0,0,0],
                      [0,10,0],
                      ],
                     [
                      [20,0,0],
                      [40,0,0]
                      ],
                     [
                      [0,0,5],
                      [0,0,15]
                      ]
                     ])
        
        expected = [30.41381265, 11.18033989,  31.6227766]
        matrix = EuclideanDistanceMatrixBuilder.calculate_geom_center(coordsets)
        numpy.testing.assert_array_almost_equal(expected, matrix.get_data())
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()