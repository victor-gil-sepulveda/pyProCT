"""
Created on 16/03/2012

@author: victor
"""
import unittest
import pyproct.clustering.test.data as test_data
from pyproct.clustering.cluster import Cluster, cluster_from_tuple, get_cluster_sizes, gen_clusters_from_class_list
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
import os


class Test(unittest.TestCase):

    def test_get_size(self):
        cluster = Cluster(prototype = 0, elements = [0,4,5,7,13])
        self.assertEqual(cluster.get_size(),5)
     
    def test_get_sizes(self):
        myclusters = []
        for c in test_data.clusters:
            myclusters.append(cluster_from_tuple(c))
        sizes = [5,4,4,4,3]
        numpy.testing.assert_array_equal(sizes, get_cluster_sizes(myclusters)[1], "Cluster sizes are different")
     
    def test_gen_clusters_from_grouping_list(self):
        #  numpy.random.random_integers(0,4,20)
        numclusters = 5
        group_list = [4, 1, 2, 2, 4, 4, 3, 4, 2, 0, 0, 3, 3, 4, 0, 3, 1, 1, 1, 2]
        true_clusters = [Cluster(0,[0,4,5,7,13]),
                         Cluster(1,[1,16,17,18]),
                         Cluster(2,[2,3,8,19]),
                         Cluster(6,[6,11,12,15]),
                         Cluster(9,[9,10,14])]
        clusters =  gen_clusters_from_class_list(group_list)
        sorted_clusters = sorted(clusters, key=lambda c: c.prototype)
        self.assertEqual(numclusters,len(sorted_clusters))
        for i in range(numclusters):
            self.assertEqual(true_clusters[i], sorted_clusters[i])
     
    def test_clusters_are_equal(self):
        clusters = [Cluster(0,[0,4,5,7,13]),
                    Cluster(2,[2,4,5,7,13]),
                    Cluster(1,[1,16,17,18]),
                    Cluster(1,[1,16,17,18]),
                    Cluster(0,[0,4,5,7,13]),
                    Cluster(1,[1,16,15,18,19]),
                    Cluster(2,[2,3,8,19])]
         
        self.assertEqual(True,clusters[0] == clusters[0])
        self.assertEqual(False,clusters[0] == clusters[1])
        self.assertEqual(False,clusters[0] == clusters[2])
        self.assertEqual(False,clusters[0] == clusters[3])
        self.assertEqual(True, clusters[0] == clusters[4])
        self.assertEqual(False,clusters[0] == clusters[5])
        self.assertEqual(False,clusters[0] == clusters[6])
         
        self.assertEqual(True, clusters[1] == clusters[1])
        self.assertEqual(False,clusters[1] == clusters[2])
        self.assertEqual(False,clusters[1] == clusters[3])
        self.assertEqual(False,clusters[1] == clusters[4])
        self.assertEqual(False, clusters[1] == clusters[5])
        self.assertEqual(False,clusters[1] == clusters[6])
         
        self.assertEqual(True,clusters[2] == clusters[2])
        self.assertEqual(True,clusters[2] == clusters[3])
        self.assertEqual(False,clusters[2] == clusters[4])
        self.assertEqual(False, clusters[2] == clusters[5])
        self.assertEqual(False,clusters[2] == clusters[6])
         
        self.assertEqual(True,clusters[3] == clusters[3])
        self.assertEqual(False,clusters[3] == clusters[4])
        self.assertEqual(False, clusters[3] == clusters[5])
        self.assertEqual(False,clusters[3] == clusters[6])
         
        self.assertEqual(True,clusters[4] == clusters[4])
        self.assertEqual(False, clusters[4] == clusters[5])
        self.assertEqual(False,clusters[4] == clusters[6])
 
         
    def test_cluster_from_tuple(self):
        cluster_tuple = (1,[16,17,18])
        expected_cluster = Cluster(1,[1,16,17,18])
        cluster = cluster_from_tuple(cluster_tuple)
        self.assertEquals(expected_cluster,cluster)
         
    def test_creation(self):
        try:
            Cluster(1,[16,17,18])
            self.fail()
        except Exception:
            pass
         
        cluster = Cluster(1,[16,1,17,18])
        cluster_copy = Cluster(1,[16,1,17,18])
        elements = [16,1,17,18]
        obtained = cluster.all_elements
         
        # Are they equal?
        numpy.testing.assert_array_equal(elements, obtained)
         
        # A modification in this list modifies the cluster
        obtained[2] = -1
        self.assertNotEquals(cluster,cluster_copy)
     
    def test_calculate_biased_medoid(self):
        condensed_matrix = CondensedMatrix([1.0, 4.5, 7.2, 6.7, 
                                                 8.5, 4.5, 3.6, 
                                                      7.8, 2.2, 
                                                           2.0]) 
        c = Cluster(None,[0,2,3,4])
        interesting_elements = [3,4,0]
        self.assertEquals(4, c.calculate_biased_medoid(condensed_matrix,interesting_elements))
        interesting_elements = [4,2,3]
        self.assertEquals(4,c.calculate_biased_medoid(condensed_matrix,interesting_elements))
     
    def test_calculate_biased_medoid_scenario(self):
        cluster = Cluster.from_dic({
                                    "prototype": 28,
                                    "elements": "0:46, 49, 51, 53, 57:58, 62:67",
                                    "id": "cluster_0"
        })
        matrix = CondensedMatrix(list(numpy.asfarray(numpy.load(os.path.join(test_data.__path__[0],"matrix.npy")))))
        self.assertEqual(cluster.prototype, cluster.calculate_medoid(matrix))
        
        cluster = Cluster.from_dic({
                                "prototype": 54,
                                "elements": "0:117, 119:135, 138:139, 141, 143, 145:146, 148:150, 153, 155:156, 167:168, 170:172, 175, 177, 190, 193, 212, 215, 234",
                                "id": "cluster_0"
                            })
        self.assertEqual(cluster.prototype, cluster.calculate_medoid(matrix))
        
        cluster = Cluster.from_dic({
                        "prototype": 1604,
                        "elements": "224, 290, 312, 334, 378, 422, 444, 466, 468, 488, 504, 526, 645, 782, 799, 821, 843, 953, 1208, 1254, 1276, 1291, 1313, 1320, 1357, 1445, 1450, 1467, 1472, 1489, 1494, 1516, 1538, 1560, 1582, 1591, 1604, 1613, 1626, 1635, 1671, 1693, 1767, 1789, 1811, 1833, 1841, 1855, 1877, 1899, 1921, 1943, 1965, 2007, 2049, 2070, 2091, 2112, 2203",
                        "id": "cluster_18"
                    })
        self.assertEqual(cluster.prototype, cluster.calculate_medoid(matrix))
        
    def test_random_sample(self):
        cluster = Cluster(None, range(0,100))
         
        self.assertItemsEqual(cluster.get_random_sample(10, 123), [45, 66, 89, 62, 67, 51, 65, 56, 22, 77])
         
    def test_to_dic(self):
        true_clusters = [Cluster(0,[0,4,5,7,13]),
                         Cluster(1,[1,16,17,18]),
                         Cluster(2,[2,3,8,19]),
                         Cluster(6,[6,11,12,15]),
                         Cluster(9,[9,10,14])]
         
        dic_clusters = [
                        {'prototype': 0, 'elements': '0, 4:5, 7, 13'},
                        {'prototype': 1, 'elements': '1, 16:18'},
                        {'prototype': 2, 'elements': '2:3, 8, 19'},
                        {'prototype': 6, 'elements': '6, 11:12, 15'},
                        {'prototype': 9, 'elements': '9:10, 14'}
                        ]
         
        for i in range(len(true_clusters)):
            self.assertDictEqual(Cluster.to_dic(true_clusters[i]), dic_clusters[i])
             
    def test_from_dic(self):
        clusters = [
                    {
                        "prototype": 400,
                        "elements": "400:410, 0, 1 ,2,3"
                    },
                    {
                        "prototype": 500,
                        "elements": "4,500:510, 5, 6:10, 11"
                    }
                ]
         
        expected_elements =[
                            [400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 0, 1, 2, 3],
                            [4, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 5, 6, 7, 8, 9, 10, 11]
                            ]
         
        for i in range(len(clusters)):
            self.assertEqual(Cluster.from_dic(clusters[i]).all_elements, expected_elements[i])
                    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()