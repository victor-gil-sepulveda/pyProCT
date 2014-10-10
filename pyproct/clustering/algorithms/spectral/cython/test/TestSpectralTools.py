"""
Created on Oct 10, 2014

@author: victor
"""
import unittest

import pyproct.clustering.algorithms.spectral.cython.spectralTools as SpectralTools
from pyRMSD.condensedMatrix import CondensedMatrix

class TestSpectralTools(unittest.TestCase):

    def test_calculate_degree_matrix(self):
        """
        Test provided by Nancy-Sarah Yacovzada.
        """
        
#         0--3
#         |\  \ 
#         5 4  2
#         
#         X = np.array(
#              [[0.0, 0.0, 1.0, 1.0, 1.0],
#               [0.0, 0.0, 1.0, 0.0, 0.0],
#               [1.0, 1.0, 0.0, 0.0, 0.0],
#               [1.0, 0.0, 0.0, 0.0, 0.0],
#               [1.0, 0.0, 0.0, 0.0, 0.0]]
#               )
        data = [0.0, 1.0, 1.0, 1.0,
                     1.0, 0.0, 0.0,
                          0.0, 0.0,
                               0.0]
        
        W = CondensedMatrix(data)
        
        self.assertSequenceEqual([ 3., 1., 2., 1., 1.], SpectralTools.calculate_degree_matrix(W))
        
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_degree']
    unittest.main()