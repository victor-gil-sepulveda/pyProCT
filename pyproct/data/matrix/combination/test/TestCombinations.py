'''
Created on 27/11/2014

@author: victor
'''
import numpy
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.data.matrix.combination.combinationMatrixCalculator import combine

class TestCombinations(unittest.TestCase):

    def test_do_combine(self):
        a = numpy.array([ 1., 4.,
                              7.])
        
        b = numpy.array([ 2., 3., 
                              1.])

        matrices = {
                    "a": CondensedMatrix(a),
                    "b": CondensedMatrix(b)
        }
        
        # Test sum
        operation = ["add","a","b"]
        c = combine(operation, matrices)
        numpy.testing.assert_array_equal( c.get_data(), [ 3.,  7.,  8.])
        
        # Test subtraction
        operation = ["sub","a","b"]
        c = combine(operation, matrices)
        numpy.testing.assert_array_equal( c.get_data(), [-1.,  1.,  6.])
        
        # Test scalar mult
        operation = ["mult",2,"a"]
        c = combine(operation, matrices)
        numpy.testing.assert_array_equal( c.get_data(), [  2.,   8.,  14.])
        
        # Test scalar mult conmutativity
        operation = ["mult","a",2]
        c = combine(operation, matrices)
        numpy.testing.assert_array_equal( c.get_data(), [  2.,   8.,  14.])
        
        # Test error conditions
        operation = ["mult",2,2]
        print "Error means OK-> ",
        self.assertRaises(SystemExit, combine, operation, matrices)
        operation = ["mult","a","a"]
        print "Error means OK-> ",
        self.assertRaises(SystemExit, combine, operation, matrices)
        operation = ["mult","dummy","a"]
        print "Error means OK-> ",
        self.assertRaises(SystemExit, combine, operation, matrices)

        # Test recursivity
        operation = ["mult",["add","a","b"], 2]
        c = combine(operation, matrices)
        numpy.testing.assert_array_equal( c.get_data(),[  6., 14.,  16.] )

        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()