'''
Created on 21/03/2012

@author: victor
'''
import unittest
from pyproclust.clustering.clusterization import Clusterization
from pyproclust.clustering.cluster import Cluster
from pyproclust.clustering.analysis.analysis import Analysis
from pyproclust.clustering.analysis.analysisRunner import AnalysisRunner



def function_without_params(clusterization):
    """
    Returns the number of clusters in a clusterization.
    """
    return str(len(clusterization.clusters))

def function_with_params(clusterization, smaller_than):
    """
    Returns how much clusters are smaller than the values in the 
    list passed as parameter.
    """
    smaller_string = ""
    for s in smaller_than:
        smaller_count = 0
        for c in clusterization.clusters:
            if c.get_size() <s:
                smaller_count += 1
        smaller_string += str(smaller_count)+" "
    return smaller_string

class Test(unittest.TestCase):

    def test_analysis_creation_and_run(self):
        clusters1 =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        
        clusters2 =(
                  Cluster(16,[16,17,18]),
                  Cluster(4,[4,5,6,7,8,15,17,18]),
                  Cluster(0,[0,1,2,3,14,15]),
                  Cluster(9,[9,10,11,12,13,14,15,14,14])
                  )
        
        clusterization1 = Clusterization(clusters1)
        clusterization2 = Clusterization(clusters2)
        
        analysis = Analysis("First analysis", function_without_params)
        analysis.run(clusterization1)
        analysis.run(clusterization2)
        self.assertEqual("First analysis\t4\t4\t",analysis.results_string)
        
        param_analysis = Analysis("Second analysis", function_with_params,[7,2,5])
        param_analysis.run(clusterization1)
        param_analysis.run(clusterization2)
        self.assertEqual("Second analysis\t3 1 2 \t2 0 1 \t",param_analysis.results_string)
        
    def test_run_and_gen_report(self):
        
        clusters1 =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        
        clusters2 =(
                  Cluster(16,[16,17,18]),
                  Cluster(4,[4,5,6,7,8,15,17,18]),
                  Cluster(0,[0,1,2,3,14,15]),
                  Cluster(9,[9,10,11,12,13,14,15,14,14])
                  )
        
        clusterization1 = Clusterization(clusters1)
        clusterization2 = Clusterization(clusters2)
        
        analyzer = AnalysisRunner()
        analysis = Analysis("First analysis", function_without_params)
        analyzer.add_analysis(analysis)
        param_analysis = Analysis("Second analysis", function_with_params,[7,2,5])
        analyzer.add_analysis(param_analysis)
        analysis = Analysis("Third analysis", function_without_params)
        analyzer.add_analysis(analysis)
        
        analyzer.run_analysis_for(clusterization1)
        analyzer.run_analysis_for(clusterization2)
        self.assertEqual("First analysis\t4\t4\t\nSecond analysis\t3 1 2 \t2 0 1 \t\nThird analysis\t4\t4\t\n", analyzer.generate_report())

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_analysis_creation']
    unittest.main()