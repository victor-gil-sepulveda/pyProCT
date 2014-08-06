"""
Created on 23/12/2013

@author: victor
"""
import unittest
from pyproct.clustering.cluster import Cluster
from pyproct.clustering.comparison.comparator import Separator, getAllElements


class Test(unittest.TestCase):

    def test_decompose(self):
        traj_ranges = {"traj_A":(0,6),"traj_B":(7,15)}
        clusters = [
                    Cluster(None,[0,1,2]),
                    Cluster(None,[3,8,10]),
                    Cluster(None,[14,4,15]),
                    Cluster(None,[5,6]),
                    Cluster(None,[7,9,11,12,13]),
                    ]

        for i in range(len(clusters)):
            clusters[i].id = str(i)

        decomposed = Separator.decompose(clusters, traj_ranges)
        all_elements = []

        for cluster_id in decomposed:
            all_elements.extend(getAllElements(decomposed[cluster_id]))

        expected = {
                    '0': {
                          'traj_A': [0, 1, 2]
                          },
                    '1': {
                          'traj_A': [3],
                          'traj_B': [8, 10]
                          },
                    '2': {
                          'traj_A': [4],
                          'traj_B': [14, 15]
                          },
                    '3': {
                          'traj_A': [5, 6]
                          },
                    '4': {
                          'traj_B': [9, 11, 12, 13, 7]
                          }
                    }
        self.assertItemsEqual(range(16),sorted(all_elements))
        self.assertDictEqual(expected, decomposed )

    def test_classify(self):
        decomposed = {
                    '0': {
                          'traj_A': [0, 1, 2]
                          },
                    '1': {
                          'traj_A': [3],
                          'traj_B': [8, 10]
                          },
                    '2': {
                          'traj_A': [4],
                          'traj_B': [14, 15]
                          },
                    '3': {
                          'traj_A': [5, 6]
                          },
                    '4': {
                          'traj_B': [9, 11, 12, 13, 7]
                          }
                    }

        expected = {
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
        self.assertDictEqual(expected, Separator.classify(decomposed))

    def test_separate(self):
        traj_ranges = {"traj_A":(0,6),"traj_B":(7,15)}
        clusters = [
                    Cluster(None,[0,1,2]),
                    Cluster(None,[3,8,10]),
                    Cluster(None,[14,4,15]),
                    Cluster(None,[5,6]),
                    Cluster(None,[7,9,11,12,13]),
                    ]

        for i in range(len(clusters)):
            clusters[i].id = str(i)

        expected = {
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
        self.assertDictEqual(expected,Separator.separate(clusters, traj_ranges))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_decompose']
    unittest.main()