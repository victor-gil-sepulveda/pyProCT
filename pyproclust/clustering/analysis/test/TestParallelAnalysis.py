'''
Created on 06/06/2012

@author: victor
'''
import unittest
from Queue import Queue
from pyproclust.clustering.analysis.parallelAnalysis import analysis_to_parallel_analysis
from pyproclust.clustering.analysis.analysis import Analysis

def useless_function(clustering):
    return 5

class Test(unittest.TestCase):

    def test_creation(self):
        analysis = Analysis("test", useless_function)
        panalysis = analysis_to_parallel_analysis(analysis)
        self.assertEqual(analysis.name, panalysis.name)
        self.assertEqual(analysis.results_string, panalysis.results_string)
        self.assertEqual(analysis.analysis_function, panalysis.analysis_function)
        self.assertEqual(analysis.number_of_results, panalysis.number_of_results)
        self.assertEqual(analysis.other_params, panalysis.other_params)

    def test_run(self):
        q = Queue()
        analysis = Analysis("test", useless_function)
        panalysis = analysis_to_parallel_analysis(analysis)
        panalysis.run("clusterization",q,10)
        self.assertEqual(1,panalysis.number_of_results)
        self.assertEqual((10,5),q.get())

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_creation']
    unittest.main()