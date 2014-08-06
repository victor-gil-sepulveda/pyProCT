"""
Created on 06/06/2012

@author: victor
"""
import unittest
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from scipy.spatial.distance import pdist

class TestKMedoids(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.condensed_matrix = CondensedMatrix([1.0, 4.5, 7.2, 6.7, 
                                                     8.5, 4.5, 3.6, 
                                                          7.8, 2.2, 
                                                               2.0]) 
        cls.kmed_alg = KMedoidsAlgorithm(cls.condensed_matrix)
        
    def test_random_seeding(self):
        for i in range(100): #@UnusedVariable
            medoids = self.kmed_alg.random_seeding(3)
            med_set = set(medoids)
            self.assertEquals(len(medoids),len(med_set))
            self.assertLessEqual(max(medoids),4)
     
    def test_convergence(self):
        list1 = [3,6,8]
        list2 = [6,3,8]
        list3 = [1,2,3]
        self.assertTrue(self.kmed_alg.convergence(list1,list2))
        self.assertFalse(self.kmed_alg.convergence(list1,list3))
     
    def test_get_closer_medoid(self):
        medoids = [2,4]
        self.assertEquals(2, self.kmed_alg.get_closer_medoid(0, medoids, self.condensed_matrix))
        self.assertEquals(4, self.kmed_alg.get_closer_medoid(3, medoids, self.condensed_matrix))
     
    def test_cluster_update(self):
        medoids = [2,4]
        self.kmed_alg.cluster_update( medoids, self.condensed_matrix)
        self.assertItemsEqual([0, 1, 0, 1, 1], self.kmed_alg.class_list)
 
    def test_gen_medoid_to_cluster_id_map(self):
        medoids = [2,4]
        self.assertEqual({2:0,4:1}, self.kmed_alg.gen_medoid_to_cluster_id_map(medoids))
     
    def test_update_medoids(self):
        points = [(0,0),(0,1),(0,-1),(1,0),
                  (6,0),(7,0),(7,1),(7,-1)]
         
        matrix = CondensedMatrix(pdist(points))
        kmed_alg = KMedoidsAlgorithm(matrix)
        kmed_alg.class_list = [0, 0, 0, 0, 
                               1, 1, 1, 1]
        numpy.testing.assert_array_equal( kmed_alg.update_medoids(), [0,5])
         
    def test_equidistant_seeding(self):
        medoids =  self.kmed_alg.equidistant_seeding(10, 100)
        expected = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95]
        numpy.testing.assert_array_equal(medoids, expected)
         
    def test_gromos_seeding(self):
        points = [(0,0),(0,1),(0,-1),(1,0),
                  (3,0),(3,1),
                  (6,0),(7,0),(7,1),(7,-1)]
#          Maximum distance of connected components is 1, for 1.5 it may discover 3 clusters.
#          For 3.2 it may discover only 2
#
#         1       5         8
#         |       |         |
#         0 - 3   4     6 - 7 
#         |                 |
#         2                 9
#         
         
        matrix = CondensedMatrix(pdist(points))
        kmed_alg = KMedoidsAlgorithm(matrix, rand_seed = 10)
         
        # With a small cutoff we get all 3 connected components
        numpy.testing.assert_array_equal( kmed_alg.gromos_seeding(3, 1.4),[0, 7, 4]) # if it's 1.5 we'll have 6 instead of 7 as medoid (as it is the first one to have 3 neighbours)
         
        # With a bigger cutoff it has to try to find the minimum cutoff for 3 clusters, then 6 is returned instead of 7
        numpy.testing.assert_array_equal( kmed_alg.gromos_seeding(3, 3.2),[3, 6, 5])
         
        # This one is regression
        numpy.testing.assert_array_equal( kmed_alg.gromos_seeding(2, 3.2),[4, 7])
         
        # This one should return a random sequence, so is only testable because of the rand_seed
        numpy.testing.assert_array_equal(kmed_alg.gromos_seeding(2, 0), [5, 3])
        
    def test_naive_case(self):
#         1       5         8
#         |       |         |
#         0 - 3   4     6 - 7 
#         |                 |
#         2                 9
        points = [(0,0),(0,1),(0,-1),(1,0),
                  (3,0),(3,1),
                  (6,0),(7,0),(7,1),(7,-1)]
        matrix = CondensedMatrix(pdist(points))
        s_algo = KMedoidsAlgorithm(matrix, 10)
        clusters = s_algo.perform_clustering({'k':3, 'seeding_type':'RANDOM'}).clusters
        
        for c in clusters:
            self.assertIn(c.prototype, [0, 4, 6])
            self.assertIn(c.all_elements, [[0, 1, 2, 3],[6, 7, 8, 9],[4, 5]])
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()