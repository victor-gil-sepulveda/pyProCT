"""
Created on 12/06/2012

@author: victor
"""
import unittest
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
import scipy.spatial.distance
from pyproct.postprocess.actions.confSpaceComparison.tools import calculate_mean_center_differences,\
    calculate_distance_stats

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

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_Name']
    unittest.main()