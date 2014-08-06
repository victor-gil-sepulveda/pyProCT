"""
Created on 23/12/2013

@author: victor
"""
import unittest
from pyproct.clustering.comparison.comparator import Analyzer
import scipy.spatial.distance
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.separated_decomposed_clusters = {
                                            "mixed":{
                                                "cluster_1": {
                                                     "traj_A":[0,1,2,3,4],
                                                     "traj_B":[5,6,7,8,9],
                                                     "traj_C":[10,11,12,13,14]
                                                }
                                            },
                                           "pure":{
                                               "cluster_2": {
                                                     "traj_A":[0,1,2,3,4]
                                               },
                                               "cluster_3": {
                                                    "traj_B":[5,6,7,8,9]
                                               }
                                            }
                                        }

        # 4 points forming a square with another point in its center
        square_points = numpy.array([[0,0], [0,2], [2,0], [2,2], [1,1]])

        # move the square to the right and up-right
        square_points_2 = square_points+numpy.array([0,5])
        square_points_3 = square_points+numpy.array([5,0])

        cls.square_points = square_points.tolist()
        cls.square_points.extend(square_points_2.tolist())
        cls.square_points.extend(square_points_3.tolist())
        cls.matrix = CondensedMatrix(scipy.spatial.distance.pdist(cls.square_points))

    def test_analyze_clustering(self):
        separated_decomposed_clusters = {
            'mixed': {
                      '1': {
                            'traj_A': [3],
                            'traj_B': [8, 10]
                            },
                      '2': {
                            'traj_A': [4],
                            'traj_B': [14, 15]
                            }
                      },
            'pure': {
                     '0': {
                           'traj_A': [0, 1, 2]
                           },
                     '3': {
                           'traj_A': [5, 6]
                           },
                     '4': {
                           'traj_B': [9, 11, 12, 13, 7]
                           }
                     }
        }

        expected_analysis = {
                             'num_pure': 3,
                             'num_mixed': 2,
                             'num_mixed_elements': 6,
                             'total_num_clusters': 5,
                             'num_pure_elements': 10,
                             'total_num_elements': 16
                             }
        analysis = {}
        Analyzer.analyze_clustering(separated_decomposed_clusters, analysis)
        self.assertDictEqual(expected_analysis, analysis)

    def test_analyze_clusters(self):

        expected_analysis = {
            'cluster_2': {
                          'global': {
                                     'mean': 1.1313708305358887,
                                     'std': 0.56568541526794436,
                                     'max': 1.4142135381698608,
                                     'num_elements': 5
                                     }
                          },
            'cluster_3': {
                          'global': {
                                     'mean': 1.1313708305358887,
                                     'std': 0.56568541526794436,
                                     'max': 1.4142135381698608,
                                     'num_elements': 5
                                     }
                          },
            'cluster_1': {
                          'centers_mean_diff': 5.6903559366861982,
                          'global': {
                                     'mean': 3.3646855751673379,
                                     'std': 1.5095995219901064,
                                     'max': 5.385164737701416,
                                     'num_elements': 15
                                     },
                          'traj_A': {
                                     'mean': 1.1313708305358887,
                                     'std': 0.56568541526794436,
                                     'max': 1.4142135381698608,
                                     'num_elements': 5
                                     },
                          'traj_B': {
                                     'mean': 1.1313708305358887,
                                     'std': 0.56568541526794436,
                                     'max': 1.4142135381698608,
                                     'num_elements': 5
                                     },
                          'traj_C': {
                                     'mean': 1.1313708305358887,
                                     'std': 0.56568541526794436,
                                     'max': 1.4142135381698608,
                                     'num_elements': 5
                                     }
                          }
            }

        analysis = {}
        Analyzer.analyze_clusters(self.separated_decomposed_clusters, self.matrix, analysis)
        self.assertDictEqual(expected_analysis, analysis)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_']
    unittest.main()