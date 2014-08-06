"""
Created on 07/08/2012

@author: victor
"""
import unittest
from pyproct.clustering.filtering.clusteringFilter import ClusteringFilter

class ClusteringMock(object):
    def __init__(self, number_of_clusters, number_of_elements):
        self.number_of_clusters = number_of_clusters
        self.total_number_of_elements = number_of_elements
        self.clusters = range(number_of_clusters)
    
    def eliminate_noise(self, max_noise):
        pass
    
    def __eq__(self, other):
        return self.number_of_clusters == other.number_of_clusters and self.total_number_of_elements == other.total_number_of_elements
    
    def __neq__(self, other):
        return not self == other

class MatrixHandlerMock():
    def __init__(self,number_of_elements):
        self.distance_matrix = MatrixMock(number_of_elements)
        
class MatrixMock(object):
    def __init__(self, number_of_elements):
        self.row_length = number_of_elements
        
class TestFiltering(unittest.TestCase):
    
    def test_check_num_clusters_in_range(self):
        myFilter = ClusteringFilter({
                                        "minimum_clusters": 5,
                                        "maximum_clusters": 30,
                                     },
                                    MatrixHandlerMock(1000))
        
        self.assertItemsEqual(myFilter.check_num_clusters_in_range(ClusteringMock(number_of_clusters = 10, number_of_elements = 1000)), [])
                              
        self.assertItemsEqual( myFilter.check_num_clusters_in_range(ClusteringMock(number_of_clusters = 2, number_of_elements = 1000)),
                               [
                                {
                                 'reason': 'TOO_FEW_CLUSTERS', 
                                 'data': {'current': 2, 
                                          'minimum': 5}
                                 }
                                ]
                               )
        self.assertItemsEqual( myFilter.check_num_clusters_in_range(ClusteringMock(number_of_clusters = 35, number_of_elements = 1000)),
                               [
                                {
                                 'reason': 'TOO_MUCH_CLUSTERS', 
                                 'data': {
                                          'current': 35, 
                                          'maximum': 30
                                          }
                                 }
                                ])
        
    def test_check_noise_level(self):
        myFilter = ClusteringFilter({
                                        "maximum_noise": 15,
                                     },
                                    MatrixHandlerMock(1000))
        # 10% noise
        self.assertItemsEqual( myFilter.check_noise_level(ClusteringMock(number_of_clusters = 10, number_of_elements = 900)), [])
        
        # 15% noise
        self.assertItemsEqual( myFilter.check_noise_level(ClusteringMock(number_of_clusters = 10, number_of_elements = 850)), [])
        
        # 20% noise
        self.assertItemsEqual( myFilter.check_noise_level(ClusteringMock(number_of_clusters = 10, number_of_elements = 800)),
                               [
                                {
                                 'reason': 'TOO_MUCH_NOISE', 
                                 'data': {
                                          'current': 20.0, 
                                          'maximum': 15
                                          }
                                 }
                                ])
        
        
    def test_check_clustering(self): 
        myFilter = ClusteringFilter({
                                        "maximum_noise": 15,
                                        "maximum_clusters": 30,
                                        "minimum_clusters": 5,
                                        "minimum_cluster_size": 50
                                     },
                                    MatrixHandlerMock(1000))
        
        self.assertItemsEqual( myFilter.check_clustering(ClusteringMock(number_of_clusters = 25, number_of_elements = 900)),[])
        
        self.assertItemsEqual( myFilter.check_clustering(ClusteringMock(number_of_clusters = 50, number_of_elements = 800)),
                                [
                                 {
                                  'reason': 'TOO_MUCH_CLUSTERS', 
                                  'data': {
                                           'current': 50, 
                                           'maximum': 30
                                           }
                                  }, 
                                 {
                                  'reason': 'TOO_MUCH_NOISE', 
                                  'data': {
                                           'current': 20.0, 
                                           'maximum': 15
                                           }
                                  }
                                 ])
    
    def test_filter_repeated(self):
        clustering_info ={"clustering 1":{
                                       "clustering":ClusteringMock(number_of_clusters = 50, number_of_elements = 800)
                                       },
                          "clustering 2":{
                                       "clustering":ClusteringMock(number_of_clusters = 25, number_of_elements = 900)
                                       },
                          "clustering 3":{
                                       "clustering":ClusteringMock(number_of_clusters = 25, number_of_elements = 900)
                                       }
                          }

        myFilter = ClusteringFilter({},MatrixHandlerMock(1000))
        sel, not_sel = myFilter.filter_repeated(clustering_info,{})
        self.assertItemsEqual(sel.keys(),["clustering 1","clustering 3"])
        self.assertItemsEqual(not_sel.keys(),["clustering 2"])
        self.assertDictEqual(not_sel["clustering 2"]["reasons"][0],
                             {'reason': 'EQUAL_TO_OTHER_CLUSTERING', 'data': {'id': 'clustering 3'}})
        
    def test_filter(self):
        myFilter = ClusteringFilter({
                                        "maximum_noise": 15,
                                        "maximum_clusters": 30,
                                        "minimum_clusters": 5,
                                        "minimum_cluster_size": 50
                                     },
                                    MatrixHandlerMock(1000))
        
        clustering_info ={"clustering 1":{
                                       "clustering":ClusteringMock(number_of_clusters = 50, number_of_elements = 800)
                                       },
                          "clustering 2":{
                                       "clustering":ClusteringMock(number_of_clusters = 25, number_of_elements = 900)
                                       },
                          "clustering 3":{
                                       "clustering":ClusteringMock(number_of_clusters = 25, number_of_elements = 900)
                                       },
                          "clustering 4":{
                                       "clustering":ClusteringMock(number_of_clusters = 31, number_of_elements = 900)
                                       }
                          }
        
        selected, not_selected =  myFilter.filter(clustering_info)
        self.assert_(len(selected) == 1 and len(not_selected) == 3)
        self.assertItemsEqual(selected.keys() ,      [ "clustering 3"])
        self.assertItemsEqual(not_selected.keys() ,  ["clustering 1","clustering 2","clustering 4",])
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()