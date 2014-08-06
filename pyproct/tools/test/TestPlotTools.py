"""
Created on 13/12/2012

@author: victor
"""
import unittest
import numpy
import pyproct.tools.plotTools

class Test(unittest.TestCase):

    def test_normalize(self):
        test_list = [-3,-2,-1,0,1,2,3]
        rescaled = [-2.0, -1.333, -0.666, 0.0, 0.666, 1.333, 2.0]
        numpy.testing.assert_array_almost_equal(pyproct.tools.plotTools.normalize(test_list, 2), rescaled, 3)
   
    def test_normalize_in_range(self):
        mylist =[255,40,35,2,245]
        norm = pyproct.tools.plotTools.normalize_in_range(mylist,0.3,1)
        numpy.testing.assert_array_almost_equal(norm,[1.0, 0.409, 0.396, 0.305, 0.972],3)
    
    def test_remove_zeros(self):
        mylist = [1,0,1,0,0,0,2,2,0,2,3,4]
        expected = [1,1,2,2,2,3,4]
        numpy.testing.assert_array_almost_equal(pyproct.tools.plotTools.remove_zeros(mylist), expected)
        
    def test_tuple_to_int(self):
        self.assertItemsEqual(pyproct.tools.plotTools.tuple_to_int((0.1,0.6,7.3,-3.2)), (0,0,7,-3))
    
    def test_shorten_name(self):
        text = "En un lugar de la mancha de cuyo nombre no quiero acordarme"
        self.assertEqual(pyproct.tools.plotTools.shorten_name(text, max_length = 5), "...darme")
        self.assertEqual(pyproct.tools.plotTools.shorten_name(text), "... acordarme")
        
    def test_shrink_matrix(self):
        self.fail("TODO")
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()