'''
Created on 06/06/2012

@author: victor
'''
import unittest
from pyproclust.clustering.analysis.parallelAnalysisRunner import ParallelAnalysisRunner
import numpy
def useless_function(clustering):
    pass

class Test(unittest.TestCase):

    def test_add_analysis(self):
        self.fail("TODO")
    
    def test_run_analysis_for(self):
        self.fail("TODO")
    
    def test_queue_unpacking(self):
        self.fail("TODO")
    
    def test_numerical_results_normalization(self):
        arr = [1,2,3,4,5,6,7,8,9,10]
        expected_result = [ 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1. ]
        prunner = ParallelAnalysisRunner(100)
        numpy.testing.assert_array_equal(expected_result, prunner.numerical_results_normalization(arr))
    
    def test_generate_report(self):
        self.fail("TODO")
    
    def test_serial_gen_report(self):
        self.fail("TODO")
        
    def test_serial_result_packing(self):
        self.fail("TODO")
        
    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_creation']
    unittest.main()    