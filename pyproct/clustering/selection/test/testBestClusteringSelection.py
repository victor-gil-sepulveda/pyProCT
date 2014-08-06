"""
Created on 07/09/2012

@author: victor
"""
import unittest
from pyproct.clustering.selection.bestClusteringSelector import BestClusteringSelector


class Test(unittest.TestCase):

    def test_normalize_one_evaluation(self):
        clustering_info = {
                           "clustering1":{
                                          "type":"algorithm2",
                                          "clustering":"clustering2",
                                          "parameters":{},
                                          "evaluation":{
                                                        "myeval":0.5,
                                                        "myothereval":1.0
                                                       }
                            },
                           "clustering2":{
                                          "type":"algorithm1",
                                          "clustering":"clustering1",
                                          "parameters":{},
                                          "evaluation":{
                                                        "myeval":-0.5,
                                                        "myothereval":2.0
                                                        }
                            },
                           "clustering3":{
                                          "type":"algorithm1",
                                          "clustering":"clustering3",
                                          "parameters":{},
                                          "evaluation":{
                                                        "myeval":0.2,
                                                        "myothereval":0.0
                                                        }
                            }
        }
        values_for_myeval = BestClusteringSelector.get_values_for_evaluation_type("myeval", clustering_info)
        self.assertDictEqual(values_for_myeval, {'clustering1': 0.5, 'clustering3': 0.2, 'clustering2': -0.5})
        BestClusteringSelector.normalize_one_evaluation_type("myeval", clustering_info)
        values_for_myeval = BestClusteringSelector.get_values_for_evaluation_type("Normalized_myeval", clustering_info)
        self.assertDictEqual( values_for_myeval, {'clustering1': 1.0, 'clustering3': 0.7, 'clustering2': 0.0})
    
    def test_get_score_for_criteria(self):
        criteria = {
                    "analysis_1":{
                                  "action": ">",
                                  "weight": 1.0
                                  },
                    "analysis_2":{
                                  "action": "<",
                                  "weight": 0.5
                                  }
        }
        
        clustering_info = {
                           'Clustering 1': {
                                             'evaluation': {
                                                            'Normalized_analysis_1': 1,
                                                            'Normalized_analysis_2': 0.3
                                              }   
                            },
                            'Clustering 2': {
                                             'evaluation': {
                                                            'Normalized_analysis_1': 0.7,
                                                            'Normalized_analysis_2': 1.
                                              }
                            }, 
                            'Clustering 3': {
                                             'evaluation': {
                                                            'Normalized_analysis_1': 0.6,
                                                            'Normalized_analysis_2': 0.5
                                              }
                            }, 
                            'Clustering 4': {
                                             'evaluation': {
                                                            'Normalized_analysis_1': 0.9,
                                                            'Normalized_analysis_2': 0.0
                                              }
                            } 
        }
        self.assertEqual(BestClusteringSelector.get_score_for_criteria("Clustering 1", 
                                                                      clustering_info,
                                                                      criteria), 
                         0.9)
        
        self.assertAlmostEqual(BestClusteringSelector.get_score_for_criteria("Clustering 4", 
                                                                      clustering_info,
                                                                      criteria), 
                         0.933333333333,12)
        
    def test_get_scores_for_all_clusters_and_criterias(self):
        criteria = {
                     "criteria 1":{
                                "analysis_1":{
                                              "action": ">",
                                              "weight": 1.0
                                              },
                                "analysis_2":{
                                              "action": "<",
                                              "weight": 0.5
                                              }
                    },
                    "criteria 2":{
                                "analysis_1":{
                                              "action": ">",
                                              "weight": 0.4
                                              },
                                "analysis_2":{
                                              "action": "<",
                                              "weight": 0.2
                                              }
                    }
        }
        
        clustering_info = {
                           'Clustering 1': {
                                             'evaluation': {
                                                            'Normalized_analysis_1': 1,
                                                            'Normalized_analysis_2': 0.3
                                              }   
                            },
                            'Clustering 2': {
                                             'evaluation': {
                                                            'Normalized_analysis_1': 0.7,
                                                            'Normalized_analysis_2': 1.
                                              }
                            }, 
                            'Clustering 3': {
                                             'evaluation': {
                                                            'Normalized_analysis_1': 0.6,
                                                            'Normalized_analysis_2': 0.5
                                              }
                            }, 
                            'Clustering 4': {
                                             'evaluation': {
                                                            'Normalized_analysis_1': 0.9,
                                                            'Normalized_analysis_2': 0.0
                                              }
                            } 
        }
        
        # regression, checked
        self.assertDictEqual({
                             'criteria 1': {
                                            'Clustering 4': 0.9333333333333332, 
                                            'Clustering 2': 0.4666666666666666, 
                                            'Clustering 3': 0.5666666666666667, 
                                            'Clustering 1': 0.9
                                            }, 
                             'criteria 2': {
                                            'Clustering 4': 0.9333333333333332, 
                                            'Clustering 2': 0.46666666666666656, 
                                            'Clustering 3': 0.5666666666666665, 
                                            'Clustering 1': 0.8999999999999999
                                            }
                             },
                             BestClusteringSelector.get_scores_for_all_clusters_and_criterias(criteria, clustering_info))
        
    def test_get_best_clustering(self):
        scores = {
                  'criteria 1': {
                                 'Clustering 4': 1.4, 
                                 'Clustering 2': 0.7, 
                                 'Clustering 3': 0.85, 
                                 'Clustering 1': 1.35}, 
                  'criteria 2': {
                                 'Clustering 4': 0.56, 
                                 'Clustering 2': 0.28, 
                                 'Clustering 3': 0.34, 
                                 'Clustering 1': 0.54
                                 }
                  }
        bclust, bcrit, scores = BestClusteringSelector.get_best_clustering(scores)
        self.assertItemsEqual( (bclust, bcrit, scores[bcrit][bclust]),  ('Clustering 4', 'criteria 1', 1.4))
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()