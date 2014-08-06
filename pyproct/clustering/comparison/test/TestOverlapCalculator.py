'''
Created on 27/02/2014

@author: victor
'''
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.clustering.comparison.overlapCalculator import OverlapCalculator
import numpy


class Test(unittest.TestCase):


    def test_calculate_global_overlap(self):
        distance_matrix = CondensedMatrix([1., 0.7, 2.,
                                               0.3, 1.,
                                                    0.7])

        decomposed_clusters = [{"traj_0":[0],"traj_1":[1]},{"traj_0":[2],"traj_1":[3]}]

        self.assertEqual(0., OverlapCalculator.calculate_global_overlap(decomposed_clusters, distance_matrix, 1, 1))

        decomposed_clusters = [{"traj_0":[0],"traj_1":[1]}, {"traj_0":[2]}, {"traj_1":[3]}]

        self.assertEqual(0., OverlapCalculator.calculate_global_overlap(decomposed_clusters, distance_matrix, 1, 1))

        #TODO: needs at least one more complex example

    def test_calculate_cluster_overlap(self):
        distance_matrix = CondensedMatrix([1., 0.7,
                                               0.3])

        decomposed_cluster = {"traj_0":[0],"traj_1":[1],"traj_2":[2]}


        self.assertAlmostEqual(0.481481488022, OverlapCalculator.calculate_cluster_overlap(1, decomposed_cluster, distance_matrix),12)
        self.assertAlmostEqual(0.4761904843, OverlapCalculator.calculate_cluster_overlap(2, decomposed_cluster, distance_matrix), 12)

        decomposed_cluster = {"traj_0":[0],"traj_1":[1]}
        self.assertAlmostEqual(1., OverlapCalculator.calculate_cluster_overlap(1, decomposed_cluster, distance_matrix), 12)
        self.assertAlmostEqual(1., OverlapCalculator.calculate_cluster_overlap(2, decomposed_cluster, distance_matrix), 12)

    def test_get_cluster_min_max_distances(self):
        # First test
        distance_matrix = CondensedMatrix([1., 0.7,
                                               0.3])

        decomposed_cluster = {"traj_0":[0],"traj_1":[1],"traj_2":[2]}

        expected_min_d, expected_max_d = [ 0.7, 0.3, 0.3], [ 1., 1., 0.7]

        min_d, max_d = OverlapCalculator.get_cluster_min_max_distances(decomposed_cluster, distance_matrix)

        numpy.testing.assert_array_almost_equal(min_d, expected_min_d, 8)
        numpy.testing.assert_array_almost_equal(max_d, expected_max_d, 8)

        #Second test
        distance_matrix = CondensedMatrix([1., 0.7, 2.,
                                               0.3, 1.,
                                                    0.7])


        decomposed_cluster = {"traj_0":[0,3],"traj_1":[1],"traj_2":[2]}

        expected_min_d, expected_max_d = [ 0.7,  0.7,  0.3,  0.3], [ 1., 1., 1.,  0.7]

        min_d, max_d =  OverlapCalculator.get_cluster_min_max_distances(decomposed_cluster, distance_matrix)

        numpy.testing.assert_array_almost_equal(min_d, expected_min_d, 8)
        numpy.testing.assert_array_almost_equal(max_d, expected_max_d, 8)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_']
    unittest.main()