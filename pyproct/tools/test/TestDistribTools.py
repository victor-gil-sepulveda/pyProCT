"""
Created on 14/05/2012

@author: victor
"""
import unittest
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.tools.distributionTools import get_distances_for_elems

class TestDistributionTools(unittest.TestCase):

    def test_get_distances_for_elems(self):
        condensed_matrix = CondensedMatrix( [0.2,  1.,  0.3,  .0, 
                                                  0.5,  0.6, 0.7,
                                                        0.9, 0.8,
                                                             0.4])
        numpy.testing.assert_almost_equal([1.0, 0.5,  0.9, 0.8], get_distances_for_elems([0,1,3,4],2,condensed_matrix),4)
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_calc_overlap']
    unittest.main()