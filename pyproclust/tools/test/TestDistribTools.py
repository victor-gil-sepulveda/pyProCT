'''
Created on 14/05/2012

@author: victor
'''
import unittest
from pyproclust.tools.distributionTools import calc_distribution_overlap,\
    desp_sign, get_nonzero_range, calculate_distribution_area,\
    translate_distribution, fit_translation_of_distributions
import numpy

class Test(unittest.TestCase):

    def test_calc_overlap(self):
        distr1 =            [0,0,0,0,0,1,2,3,2,1,0,0,1,2,4,2,0,0,0]
        distr2 =            [0,1,2,3,2,1,0,0,1,2,4,2,0,0,0,0,0,0,0]
        expected_overlap =  [0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0]
        self.assertItemsEqual(expected_overlap,calc_distribution_overlap(distr1,distr2,len(distr1)))
    
    def test_desp_sign(self):
        self.assertEquals(desp_sign(2,3),-1)
        self.assertEquals(desp_sign(2,2),1)
        self.assertEquals(desp_sign(2,-1),1)
        self.assertEquals(desp_sign(10,5),1)
        
    def test_get_nonzero_range(self):
        distr1 =            [0,0,0,0,0,1,2,3,2,1,0,0,1,2,4,2,0,0,0]
        distr2 =            [0,1,2,3,2,1,0,0,1,2,4,2,0,0,0,0,0,0,0]
        distr3 =            [1,0,0,0,0,0,1,0,0,0,0,0,0]
        distr4 =            [0,0,0,0,0,0,1,0,0,0,0,0,1]
        
        self.assertItemsEqual(get_nonzero_range(distr1), (5,15))
        self.assertItemsEqual(get_nonzero_range(distr2), (1,11))
        self.assertItemsEqual(get_nonzero_range(distr3), (0,6))
        self.assertItemsEqual(get_nonzero_range(distr4), (6,12))

    def test_calculate_distribution_area(self):
        distr1 =            [0,0,0,0,0,1,2,3,2,1,0,0,1,2,4,2,0,0,0]
        self.assertEqual(calculate_distribution_area(distr1, 1),numpy.sum(distr1))
    
    def test_translate_distribution(self):
        distr1 =            [0,0,0,0,0,1,2,3,2,1,0,0,1,2,4,2,0,0,0]
        distr2 =            [0,1,2,3,2,1,0,0,1,2,4,2,0,0,0,0,0,0,0]
        distr3 =            [0,0,0,0,0,0,0,0,0,0,1,2,3,2,1,0,0,1,2]
        self.assertItemsEqual(distr2, translate_distribution(distr1,-4))
        self.assertItemsEqual(distr3, translate_distribution(distr1,5))
        
    def test_fit_translation_of_distributions(self):
        distr1 =  [0,0,0,0,0,1,2,3,2,1,0,0,1,2,4,2,0,0,0]
        distr2 =  [0,1,2,3,2,1,0,0,1,2,4,2,0,0,0,0,0,0,0]
        distr3 =  [0,0,0,0,0,0,0,0,0,0,1,2,3,2,1,0,0,1,2]
        self.assertEqual( (18, 4),fit_translation_of_distributions(distr1,distr2,len(distr1),1))
        self.assertEqual( (18, -4),fit_translation_of_distributions(distr2,distr1,len(distr1),1))
        self.assertEqual( (12, -5),fit_translation_of_distributions(distr1,distr3,len(distr1),1))
        self.assertEqual( (12, 5),fit_translation_of_distributions(distr3,distr1,len(distr1),1))
        
        distr4 = [1,0,0,0,0,0,1,0,0,0,0,0,0]
        distr5 = [0,0,0,0,0,0,1,0,0,0,0,0,1]
        distr6 = [0,0,0,0,0,0,1,0,0,3,0,0,1]
        self.assertEqual((2,-6),fit_translation_of_distributions(distr4,distr5,len(distr4),1))
        self.assertEqual((2, -6),fit_translation_of_distributions(distr4,distr6,len(distr4),1))
        self.assertEqual((1, 6),fit_translation_of_distributions(distr6,distr4,len(distr4),1))
        
        distr7 = [3,2,1,0,0,0,0,0,0,0,0,0,0]
        distr8 = [0,0,0,0,0,0,0,0,0,0,1,2,3]
        distr9 = [0,0,0,0,5,4,3,2,1,0,0,0,0]
        distr10= [0,0,0,0,1,2,3,4,5,0,0,0,0]
        self.assertEqual((6, -4),fit_translation_of_distributions(distr7,distr9,len(distr7),1))
        self.assertEqual((6, 6),fit_translation_of_distributions(distr8,distr9,len(distr7),1))
        self.assertEqual((6, -6),fit_translation_of_distributions(distr7,distr10,len(distr7),1)) 
        self.assertEqual((6, 6),fit_translation_of_distributions(distr8,distr10,len(distr7),1))
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_calc_overlap']
    unittest.main()