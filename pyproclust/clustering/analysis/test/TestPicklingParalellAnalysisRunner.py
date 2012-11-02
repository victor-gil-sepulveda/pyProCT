'''
Created on 15/06/2012

@author: victor
'''
import unittest
from pyproclust.clustering.analysis.picklingParallelAnalysisRunner import PicklingParallelAnalysisRunner
import numpy
import pickle
import tempfile
import os

class FakeClustering():
    def __init__(self):
        self.details = "I'm a fake cluster" 
    def __str__(self):
        print self.name

class Test(unittest.TestCase):

    def test_run_analysis_for(self):
        prunner = PicklingParallelAnalysisRunner(100)
        prunner.run_analysis_for(FakeClustering())
        self.assertEqual(prunner.current_analysis, 1)
        self.assertEqual(len(prunner.evaluation_data), 1)
        self.assertEqual(len(prunner.process_manager.processes), 1)
        
    def test_recover_evaluation_data(self):
        prunner = PicklingParallelAnalysisRunner(100)
        results = [('Clustering 1', {'color': 'red', 'number': 1}), ('Clustering 3', {'color': 'green', 'number': 2}), ('Clustering 2', {'color': 'blue', 'number': 3}), ('Clustering 4', {'color': 'gray', 'number': 4})]
        
        for myobject in results:
            tmp_file_handler_descriptor , path = tempfile.mkstemp() 
            tmp_file_handler = os.fdopen(tmp_file_handler_descriptor,'w')
            pickle.dump(myobject, tmp_file_handler)
            tmp_file_handler.close()
            prunner.evaluation_data.append(path)
            
        ordered_clusterings = []
        all_results = {} 
        prunner.recover_evaluation_data(ordered_clusterings, all_results)
        
        self.assertItemsEqual(ordered_clusterings, ['Clustering 1', 'Clustering 3', 'Clustering 2', 'Clustering 4'])
        self.assertItemsEqual(all_results, {'color': ['red', 'green', 'blue', 'gray'], 'number': [1, 2, 3, 4]}) 
        
    def test_temporary_file_check(self):
        myobject = ["Clustering 1", "Clustering 3", "Clustering 2","Clustering 4"]
        tmp_file_handler = tempfile.TemporaryFile()
        pickle.dump(myobject, tmp_file_handler)
        tmp_file_handler.seek(0)
        self.assertEqual(myobject, pickle.load(tmp_file_handler))
        
    def test_repack_results(self):
        prunner = PicklingParallelAnalysisRunner(100)
        ordered_clusterings = ["Clustering 1", "Clustering 3", "Clustering 2","Clustering 4"]
        all_results = {"number":[1,2,3,4],"color":["red","green","blue","gray"]}
        expected_repack = [('Clustering 1', {'color': 'red', 'number': 1}), ('Clustering 3', {'color': 'green', 'number': 2}), ('Clustering 2', {'color': 'blue', 'number': 3}), ('Clustering 4', {'color': 'gray', 'number': 4})]
        self.assertItemsEqual(expected_repack, prunner.repack_results(ordered_clusterings, all_results))
         
    def test_gen_final_string_and_normalize(self):
        prunner = PicklingParallelAnalysisRunner(100)
        all_results = {"number":[1,2,3,4],"color":["red","green","blue","gray"]}
        expected_result_string = "color\tred\tgreen\tblue\tgray\t\nnumber\t1\t2\t3\t4\t\n"
        self.assertEqual(expected_result_string, prunner.gen_final_string_and_normalize(all_results))
        self.assertItemsEqual(all_results["number"],[ 0.25, 0.5, 0.75, 1.])
        self.assertItemsEqual(all_results["color"],["red","green","blue","gray"])

    def test_gen_results_string(self):
        prunner = PicklingParallelAnalysisRunner(100)
        self.assertEqual("number\t1\t2\t3\t4\t", prunner.gen_results_string( "number", [1,2,3,4]))
        self.assertEqual("color\tred\tgreen\tblue\t", prunner.gen_results_string( "color", ["red","green","blue"]))
    
    def test_numerical_results_normalization(self):
        arr = [1,2,3,4,5,6,7,8,9,10]
        expected_result = [ 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1. ]
        prunner = PicklingParallelAnalysisRunner(100)
        numpy.testing.assert_array_equal(expected_result, prunner.numerical_results_normalization(arr))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_generate_report']
    unittest.main()