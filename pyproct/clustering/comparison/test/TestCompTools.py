'''
Created on 03/03/2014

@author: victor
'''
import unittest
from pyproct.clustering.comparison.tools import getAllElements,\
    mergeSeparatedClusters, calculate_mean_center_differences,\
    calculate_distance_stats
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
import scipy.spatial.distance


class Test(unittest.TestCase):


    @classmethod
    def setUpClass(cls):
        cls.decomposed_cluster = {
                              "traj_A":[0,1,2,3,4],
                              "traj_B":[5,6,7,8,9],
                              "traj_C":[10,11,12,13,14]
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

    def test_calculate_mean_centers_difference(self):
        expected_medoids = [[1,1],[1,6],[6,1]]
        expected_mean = numpy.mean(scipy.spatial.distance.pdist(expected_medoids))
        self.assertAlmostEqual(expected_mean, calculate_mean_center_differences(self.decomposed_cluster, self.matrix), 8)

    def test_calculate_distance_stats(self):
        calc_mean, calc_std, calc_radius = calculate_distance_stats(self.decomposed_cluster["traj_A"], self.matrix)
        distances = [1.4142135381698608, 1.4142135381698608, 1.4142135381698608, 1.4142135381698608, 0.0]
        expected_mean = numpy.mean(distances)
        expected_std = numpy.std(distances)
        expected_radius = numpy.max(distances)
        self.assertItemsEqual((expected_mean,expected_std,expected_radius),(calc_mean, calc_std, calc_radius))

    def test_getAllElements(self):
        numpy.testing.assert_array_equal( sorted(getAllElements(self.decomposed_cluster)), range(15))

    def test_mergeSeparatedClusters(self):
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

        expected = [{
                        'traj_A': [3],
                        'traj_B': [8, 10]
                    },
                   {
                        'traj_A': [4],
                        'traj_B': [14, 15]
                    },
                    {
                       'traj_A': [0, 1, 2]
                    },
                    {
                       'traj_A': [5, 6]
                    },
                    {
                       'traj_B': [9, 11, 12, 13, 7]
                    }
                  ]

        self.assertItemsEqual(expected, mergeSeparatedClusters(separated_decomposed_clusters))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()