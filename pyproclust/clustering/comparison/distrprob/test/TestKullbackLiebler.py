'''
Created on 12/12/2012

@author: victor
'''
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproclust.clustering.comparison.distrprob.kullbackLieblerDivergence import KullbackLeiblerDivergence,\
    smoothed
import numpy

class Test(unittest.TestCase):


    def test_get_matrix_data(self):
        condensed_matrix = CondensedMatrix([1.0, 4.5, 7.2, 3.3, 6.8, 6.1, 
                                                 8.5, 4.5, 4.6, 9.0, 1.0,
                                                      7.8, 1.0, 0.0, 6.5,
                                                           9.6, 2.9, 2.2,
                                                                4.4, 7.1,
                                                                     8.0])
        
        numpy.testing.assert_almost_equal(KullbackLeiblerDivergence.get_matrix_data(condensed_matrix, 0, 5),\
                                          [ 1., 4.5, 7.2, 3.3, 8.5, 4.5, 4.6, 7.8, 1., 9.6],\
                                          5)
         
        numpy.testing.assert_almost_equal(KullbackLeiblerDivergence.get_matrix_data(condensed_matrix, 4, 3),\
                                          [4.4, 7.1,  8.],\
                                          5)
        
    def test_kullback_leibler_divergence_calculation(self):
        """
        Regression Test.
        """
        first_distribution_probs = [ 0.83, 0.1, 0.07]
        second_distribution_props = [0.65, 0.2, 0.15]
        screwed_distribution_props = [0.3, 0.0, 0.7]
        
        self.assertAlmostEqual(KullbackLeiblerDivergence.kullback_leibler_divergence_calculation(first_distribution_probs, second_distribution_props),\
                               0.115749946056,
                               10)
        
    def test_smoothing(self):
        """
        Regression Test.
        """
        screwed_distribution_props = [0.3, 0.0, 0.7]
        numpy.testing.assert_almost_equal([  2.99999995e-01,   1.00000000e-08,   6.99999995e-01,],smoothed(screwed_distribution_props),8)
        
    def test_kullback_leibler_creation_calculation(self):
        """
        Regression Test.
        """
        condensed_matrix = CondensedMatrix([1.0, 4.5, 7.2, 3.3, 6.8, 6.1, 
                                                 8.5, 4.5, 4.6, 9.0, 1.0,
                                                      7.8, 1.0, 0.0, 6.5,
                                                           9.6, 2.9, 2.2,
                                                                4.4, 7.1,
                                                                     8.0])
        kl = KullbackLeiblerDivergence("pdb1", "pdb2", 5, 3, condensed_matrix)
        numpy.testing.assert_almost_equal( kl.get_calculated_KL_values(), (23.653473577868997, 25.657101457909057), 10)
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_get_matrix_data']
    unittest.main()