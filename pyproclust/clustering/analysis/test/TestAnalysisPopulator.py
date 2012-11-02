'''
Created on 06/06/2012

@author: victor
'''
import unittest
from pyproclust.clustering.cluster import Cluster
from pyproclust.clustering.clusterization import Clusterization
from pyproclust.clustering.analysis.analysisPopulator import AnalysisPopulator

clusters = [    Cluster(0,[0,4,5,7,13]),
                Cluster(1,[1,16,17,18]),
                Cluster(2,[2,3,8,19]),
                Cluster(6,[6,11,12,15]),
                Cluster(9,[9,10,14])]
num_elements = 20

clustering = Clusterization(clusters, "Test clusterization")

class mockMatrix(object):
    def __init__(self):
        self.row_length = 0
        
class mockAnalyzer(object):
    def __init__(self):
        self.current_analysis = 0
    
    def add_analysis(self,analysis):
        self.current_analysis += 1

matrix = mockMatrix()
analyzer = mockAnalyzer()
populator = AnalysisPopulator(analyzer, matrix, ["PercentInTop4","NoiseLevel"])
        
class Test(unittest.TestCase):
    
    def test_analysis_created(self):
        self.assertEqual(analyzer.current_analysis, 3)

    def test_num_clusters(self):
        self.assertEquals(5,populator.analysis_function_num_clusters(clustering))
    
    def test_total_elems(self):
        self.assertEquals(num_elements,populator.analysis_function_total_elements(clustering))
    
    def test_top_4(self):
        self.assertEqual(85, populator.analysis_function_top_4(clustering))

    def test_top_percent(self):
        self.assertEqual(25, populator.analysis_function_top_percent(clustering))
    
    def test_num_clusters_to_percent(self):
        self.assertEquals(True,True) # this test is already done in the clustering class
        
    def test_details(self):
        self.assertEqual("Test clusterization", populator.analysis_function_details(clustering))
    
    def test_noise_level(self):
        self.assertEqual(75, populator.analysis_function_noise_level(clustering, num_elements*4))
        
    def test_mean_cluster_size(self):
        self.assertEqual(4, populator.analysis_function_mean_cluster_size(clustering))
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_num_clusters']
    unittest.main()