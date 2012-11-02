'''
Created on 06/06/2012

@author: victor
'''
import unittest
from pyproclust.matrix.condensedMatrix import CondensedDistanceMatrix
from pyproclust.algorithms.kmedoids.kMedoids import KMedoids

condensed_matrix = CondensedDistanceMatrix([1.0, 4.5, 7.2, 6.7, 
                                                 8.5, 4.5, 3.6, 
                                                      7.8, 2.2, 
                                                           2.0]) 
kmed_alg = KMedoids(condensed_matrix)

class Test(unittest.TestCase):

    def test_random_seeding(self):
        for i in range(100): #@UnusedVariable
            medoids = kmed_alg.random_seeding(3)
            med_set = set(medoids)
            self.assertEquals(len(medoids),len(med_set))
            self.assertLessEqual(max(medoids),4)
    
    def test_convergence(self):
        list1 = [3,6,8]
        list2 = [6,3,8]
        list3 = [1,2,3]
        self.assertTrue(kmed_alg.convergence(list1,list2))
        self.assertFalse(kmed_alg.convergence(list1,list3))
    
    def test_get_closer_medoid(self):
        medoids = [2,4]
        self.assertEquals(2, kmed_alg.get_closer_medoid(0, medoids, condensed_matrix))
        self.assertEquals(4, kmed_alg.get_closer_medoid(3, medoids, condensed_matrix))
    
    def test_cluster_update(self):
        medoids = [2,4]
        kmed_alg.cluster_update( medoids, condensed_matrix)
        self.assertItemsEqual([0, 1, 0, 1, 1],kmed_alg.class_list)

    def test_gen_medoid_to_cluster_id_map(self):
        medoids = [2,4]
        self.assertEqual({2:0,4:1}, kmed_alg.gen_medoid_to_cluster_id_map(medoids))
    
    def test_update_medoids(self):
        kmed_alg.class_list = [0, 1, 0, 1, 1]
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()