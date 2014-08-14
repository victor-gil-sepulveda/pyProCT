"""
Created on 06/06/2012

@author: victor
"""
import unittest
from pyproct.clustering.cluster import Cluster
from pyproct.clustering.clustering import Clustering
from pyproct.clustering.evaluation.analysis.analysisPopulator import AnalysisPopulator
        
class AnalysisPopulatorMock(AnalysisPopulator):
    def __init__(self ,parameters):
        AnalysisPopulator.__init__(self, "", "", parameters)
    
    def build_all_analysis(self, one, two):
        self.all_possible_analysis = {
                                        "Analysis1": "Analysis Object 1",
                                        "Analysis2": "Analysis Object 2",
                                        "Analysis3": "Analysis Object 3",
                                        "Analysis4": "Analysis Object 4"
        }
        
class ClusteringMock(object):
    def __init__(self, number_of_clusters, number_of_elements):
        self.number_of_clusters = number_of_clusters
        self.total_number_of_elements = number_of_elements
        self.clusters = range(number_of_clusters)
        self.details = "ClusteringMock"
    
    def sort_clusters_by_size(self):
        pass
    
    def get_population_percent_of_cluster(self, i):
        return self.clusters[i]
    
    def get_population_percent_of_n_bigger_clusters(self, n):
        return [60,20,10,5]
    
   
class TestAnalysisPopulator(unittest.TestCase):
    
    def test_get_query_and_evaluation_analysis_types(self):
        parameters = {
                      "clustering":{
                            "evaluation": {
                                            "evaluation_criteria": {
                                                                    "criteria_0": {
                                                                                   "CythonMirrorCohesion":{
                                                                                        "action": ">",
                                                                                        "weight": 0.05
                                                                                    },
                                                                                    "CythonMinimumMeanSeparation":{
                                                                                        "action": ">",
                                                                                        "weight": 0.1
                                                                                    },
                                                                                    "CythonSilhouette":{
                                                                                        "action": ">",
                                                                                        "weight": 0.15
                                                                                    }
                                                                    }
                                            },
                                            "query_types": [
                                                            "NumClusters",
                                                            "CythonMinimumMeanSeparation",
                                                            "NoiseLevel"
                                            ]
                                           }
                            }
                      }
        
        self.assertItemsEqual( AnalysisPopulator.get_evaluation_analysis_types(parameters),
            ['CythonMinimumMeanSeparation', 'CythonMirrorCohesion', 'CythonSilhouette'])

        self.assertItemsEqual( AnalysisPopulator.get_query_and_evaluation_analysis_types(parameters),
            ['CythonMinimumMeanSeparation', 'NumClusters',  'CythonMirrorCohesion', 'NoiseLevel', 'CythonSilhouette'])

    def test_populate_analyzer(self):
        queue = AnalysisPopulatorMock({
            "clustering":{
                          "evaluation": {
                            "evaluation_criteria": {
                                                    "criteria_0": {
                                                                    "Analysis1":{
                                                                        "action": ">",
                                                                        "weight": 0.05
                                                                    },
                                                                    "Analysis2":{
                                                                        "action": ">",
                                                                        "weight": 0.1
                                                                    },
                                                                    "Analysis3":{
                                                                        "action": ">",
                                                                        "weight": 0.15
                                                                    }
                                                    }
                            },
                            "query_types": [
                                            "Analysis2",
                                            "Analysis3",
                                            "Analysis4"
                            ]
                           }
            }
        }).get_analysis_list()
        
        self.assertItemsEqual(queue, ["Analysis Object 1", "Analysis Object 2", "Analysis Object 3", "Analysis Object 4"])
    
    def test_num_clusters(self):
        analysisPopulator = AnalysisPopulatorMock("")
        self.assertEquals(5,analysisPopulator.analysis_function_num_clusters(ClusteringMock(5, 1000)))
     
    def test_total_elems(self):
        analysisPopulator = AnalysisPopulatorMock("")
        self.assertEquals(1000,analysisPopulator.analysis_function_total_elements(ClusteringMock(5, 1000)))
     
    def test_top_4(self):
        analysisPopulator = AnalysisPopulatorMock("")
        self.assertEqual(95, analysisPopulator.analysis_function_top_4(ClusteringMock(5, 1000)))
 
    def test_top_percent(self):
        analysisPopulator = AnalysisPopulatorMock("")
        self.assertEqual(0, analysisPopulator.analysis_function_top_percent(ClusteringMock(5, 1000)))
      
    def test_details(self):
        analysisPopulator = AnalysisPopulatorMock("")
        self.assertEqual("ClusteringMock", analysisPopulator.analysis_function_details(ClusteringMock(5, 1000)))
     
    def test_noise_level(self):
        analysisPopulator = AnalysisPopulatorMock("")
        self.assertEqual(25, analysisPopulator.analysis_function_noise_level(ClusteringMock(5, 750),1000))
         
    def test_mean_cluster_size(self):
        clusters = [    Cluster(0,[0,4,5,7,13]),
                        Cluster(1,[1,16,17,18]),
                        Cluster(2,[2,3,8,19]),
                        Cluster(6,[6,11,12,15]),
                        Cluster(9,[9,10,14])]
        clustering = Clustering(clusters, "Test Clustering")
        analysisPopulator = AnalysisPopulatorMock("")
        self.assertEqual(4, analysisPopulator.analysis_function_mean_cluster_size(clustering))
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_num_clusters']
    unittest.main()