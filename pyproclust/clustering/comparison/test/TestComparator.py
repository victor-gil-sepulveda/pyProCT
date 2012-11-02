'''
Created on 12/06/2012

@author: victor
'''
import unittest
from pyproclust.clustering.cluster import Cluster
import scipy.spatial.distance
from pyproclust.clustering.clusterization import Clustering
from pyproclust.clustering.comparison.comparator import Separator,\
    getTotalNumberOfElements, dic2ListWoZeros, calculate_centers_difference,\
    calculate_cluster_mean_distance_and_dispersion, ClusterStatisticalData,\
    ClusteringStatisticalAnalyzer, getAllOfAnAttribute, calculate_max_distance
from pyproclust.matrix.condensedMatrix import CondensedDistanceMatrix
import math

class ClusteringStatisticalAnalyzerTest(ClusteringStatisticalAnalyzer):
    def __init__(self, clustering, pure_A, pure_B, mixed, traj_A_pdb, traj_B_pdb, condensed_matrix,traj_1_numbr_of_models,traj_2_numbr_of_models,matrix):
        self.matrix = condensed_matrix
        self.clustering = clustering
        self.total_number_of_elements = clustering.total_number_of_elements
        self.statistics_dic = {}
        self.pure_A = pure_A 
        self.pure_B = pure_B
        self.mixed_clusters_with_elements = mixed
        self.traj_A_number_of_models = traj_1_numbr_of_models
        self.traj_B_number_of_models = traj_2_numbr_of_models
        self.mixed_clusters_without_elements = []
        for m in self.mixed_clusters_with_elements:
            self.mixed_clusters_without_elements.append(m[0])
        self.matrix = matrix
        self.cluster_statistical_data =[]

class StupidClass(object):
    def __init__(self,attr1,attr2,attr3):
        self.one = attr1
        self.two = attr2
        self.three = attr3

class Test(unittest.TestCase):

    def test_separate_pure_from_mixed_clusters(self):
        traj_1_elems = [0,1,2,3,4,5,6]
        traj_2_elems = [7,8,9,10,11]
        clusters = [Cluster(None,elements = [3,4,5]),Cluster(None,elements = [8,9]),
                    Cluster(None,elements = [0,1,2,7]),Cluster(None,elements = [6,10,11])]
        clustering = Clustering(clusters)
        
        separator = Separator()
        
        pure1, pure2, mixed = separator.separate_pure_from_mixed_clusters(clustering.clusters,traj_1_elems,traj_2_elems)
        
        expected_pure1 = [Cluster(None,elements = [3,4,5])]
        expected_pure2 = [Cluster(None,elements = [8,9])]
        expected_mixed = [(Cluster(None,elements = [0,1,2,7]),[0, 1, 2],[7]),
                          (Cluster(None,elements = [6,10,11]),[6],[10, 11])]
        
        self.assertEqual(expected_pure1, pure1)
        self.assertEqual(expected_pure2, pure2)
        self.assertItemsEqual(expected_mixed, mixed)

    def test_get_total_number_of_elements(self):
        pure1 = [Cluster(None,elements = [3,4,5])]
        pure2 = [Cluster(None,elements = [8,9])]
        mixed = [Cluster(None,elements = [0,1,2,7]), Cluster(None,elements = [6,10,11])]
        
        num_elems = getTotalNumberOfElements([pure1,pure2,mixed])
        self.assertEqual(num_elems, 12)
        num_elems = getTotalNumberOfElements([pure1,pure2])
        self.assertEqual(num_elems, 5)
        
    def test_dic_2_list_wo_zeros(self):
        mydictionary = {"A":0,"B":1,"C":4,"D":0,"E":1,"F":8,"G":3,"H":0,"I":9}
        expected_list = [1,4,1,8,3,9]
        self.assertItemsEqual(dic2ListWoZeros(mydictionary), expected_list)

    def test_calculate_centers_difference(self):
        self.assertEqual(calculate_centers_difference('P',"c1","c2"), 0)
        # 2D points
        points = [(0,0),(0,1),(0,2),(0,3),(0,4),(1,0),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2),(2,3),(2,4),(3,0),(3,1),(3,2),(3,3),(3,4),(4,0),(4,1),(4,2),(4,3),(4,4)]
        # Global medoid
        expected_medoid = 12
        cluster = Cluster(None, range(len(points)))
        distance_data = scipy.spatial.distance.pdist(points)
        matrix = CondensedDistanceMatrix(distance_data)
        medoid = cluster.calculate_biased_medoid(matrix, cluster.all_elements)
        self.assertEqual(expected_medoid,medoid)
        self.assertEqual(calculate_centers_difference('M', cluster, matrix, [0,1,2,5,6,7,10,11,12], [12,13,14,17,18,19,22,23,24]),math.sqrt(8))
        self.assertEqual(calculate_centers_difference('M', cluster, matrix, [6], [2,3,4,7,8,9,12,13,14]),2)

    def test_calculate_cluster_mean_distance_and_dispersion(self):
        points = [(0,0),(0,1),(0,2),(0,3),(0,4),(1,0),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2),(2,3),(2,4),(3,0),(3,1),(3,2),(3,3),(3,4),(4,0),(4,1),(4,2),(4,3),(4,4)]
        distance_data = scipy.spatial.distance.pdist(points)
        matrix = CondensedDistanceMatrix(distance_data)
        cluster = Cluster(None, range(len(points)))
        self.assertEqual(calculate_cluster_mean_distance_and_dispersion(cluster, [1,5,6,7,11], matrix), (1.0, 0.0))
    
    def test_cluster_statistical_data(self):
        points = [(0,0),(0,1),(0,2),(0,3),(0,4),(1,0),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2),(2,3),(2,4),(3,0),(3,1),(3,2),(3,3),(3,4),(4,0),(4,1),(4,2),(4,3),(4,4)]
        distance_data = scipy.spatial.distance.pdist(points)
        matrix = CondensedDistanceMatrix(distance_data)
        cluster = Cluster(None, range(len(points)))
        d1 = ClusterStatisticalData(cluster,'P',matrix)
        expected = """[Cluster Stats. Data -\n\tType: 'Pure'\n\tSize: 25\n\tMean Dist. to centre: 1.952\n\tDispersion: 0.595]\n"""
        self.assertEqual(expected,str(d1))
        mixedCluster = Cluster(None,[6]+[2,3,4,7,8,9,12,13,14])
        d2 = ClusterStatisticalData(mixedCluster,'M',matrix,[6],[2,3,4,7,8,9,12,13,14])
        expected = """[Cluster Stats. Data -\n\ttype: 'Mixed'\n\tSize: 10\n\tMean Dist. to centre: 1.295\n\tDispersion: 0.317\n\tDist. between centers: 2.000\n\n\tA elements inside mixed:[Cluster Stats. Data -\n\tType: 'Pure'\n\tSize: 1\n\tMean Dist. to centre: 0.000\n\tDispersion: 0.000]\n\n\tB elements inside mixed:[Cluster Stats. Data -\n\tType: 'Pure'\n\tSize: 9\n\tMean Dist. to centre: 1.207\n\tDispersion: 0.207]\n]\n"""
        self.assertEqual(expected,str(d2))
        
    def test_calculate_separation_statistics(self):
        # A : 0, 1, 2, 3, 4, 5, 6
        # B : 7, 8, 9, 10, 11
        pure1 = [Cluster(None,elements = [3,4,5])]
        pure2 = [Cluster(None,elements = [8,9])]
        mixed = [Cluster(None,elements = [0,1,2,7]), Cluster(None,elements = [6,10,11])]
        mixed_w_elements = [(Cluster(None,elements = [0,1,2,7]), [0,1,2], [7]), (Cluster(None,elements = [6,10,11]), [6], [10,11])]
        s = ClusteringStatisticalAnalyzerTest(Clustering(pure1+pure2+mixed),pure1, pure2, mixed_w_elements, "A", "B", "condensed", 7, 5,"matrix")
        s.per_clustering_analytics()
        expected_dic = {}
        expected_dic["number_elements_pure_A"] = 3
        expected_dic["elems_percent_pure_A"]= 3*100/12.
        expected_dic["elems_self_percent_pure_A"]= 3*100 / 7.
        expected_dic["number_elements_pure_B"] = 2
        expected_dic["elems_percent_pure_B"]= 2*100 / 12.
        expected_dic["elems_self_percent_pure_B"]= 2*100 / 5.
        expected_dic["number_elements_mixed"] =  7
        expected_dic["elems_percent_mixed"] = 7*100 / 12.
        self.assertDictEqual(expected_dic, s.statistics_dic)
        
    def test_per_cluster_analytics(self):
        points = [(0,0),(0,1),(0,2),(0,3),(0,4),(1,0),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2),(2,3),(2,4),(3,0),(3,1),(3,2),(3,3),(3,4),(4,0),(4,1),(4,2),(4,3),(4,4)]
        distance_data = scipy.spatial.distance.pdist(points)
        matrix = CondensedDistanceMatrix(distance_data)
        pureA = [Cluster(None, elements = [0,1])]
        pureB = [Cluster(None, elements = [23,24])]
        mixed = [Cluster(None, elements = [3,4,9,8]),Cluster(None, elements = [15,16,20,21])]
        mixed_w_elements = [(Cluster(None, elements = [3,4,9,8]),[3,4],[9,8]),(Cluster(None, elements = [15,16,20,21]),[15,16],[20,21])]
        
        s = ClusteringStatisticalAnalyzerTest(Clustering(pureA+pureB+mixed),pureA, pureB, mixed_w_elements, "A", "B", "condensed", 7, 5,matrix)
        s.per_cluster_analytics()
        calculated = ""
        for cd in s.cluster_statistical_data:
            calculated += str(cd)+"\n"
        expected = """[Cluster Stats. Data -\n\tType: 'Pure'\n\tSize: 2\n\tMean Dist. to centre: 1.000\n\tDispersion: 0.000]\n\n[Cluster Stats. Data -\n\tType: 'Pure'\n\tSize: 2\n\tMean Dist. to centre: 1.000\n\tDispersion: 0.000]\n\n[Cluster Stats. Data -\n\ttype: 'Mixed'\n\tSize: 4\n\tMean Dist. to centre: 1.138\n\tDispersion: 0.195\n\tDist. between centers: 1.000\n\n\tA elements inside mixed:[Cluster Stats. Data -\n\tType: 'Pure'\n\tSize: 2\n\tMean Dist. to centre: 1.000\n\tDispersion: 0.000]\n\n\tB elements inside mixed:[Cluster Stats. Data -\n\tType: 'Pure'\n\tSize: 2\n\tMean Dist. to centre: 1.000\n\tDispersion: 0.000]\n]\n\n[Cluster Stats. Data -\n\ttype: 'Mixed'\n\tSize: 4\n\tMean Dist. to centre: 1.138\n\tDispersion: 0.195\n\tDist. between centers: 1.000\n\n\tA elements inside mixed:[Cluster Stats. Data -\n\tType: 'Pure'\n\tSize: 2\n\tMean Dist. to centre: 1.000\n\tDispersion: 0.000]\n\n\tB elements inside mixed:[Cluster Stats. Data -\n\tType: 'Pure'\n\tSize: 2\n\tMean Dist. to centre: 1.000\n\tDispersion: 0.000]\n]\n\n"""
        self.assertEqual(expected, calculated)
    
    def test_gather_all(self):
        stupid_objects = []
        for i in range(100):
            stupid_objects.append(StupidClass(i,100+i,200+i))
            
        return_list =  getAllOfAnAttribute(stupid_objects,"one")+getAllOfAnAttribute(stupid_objects,"two")+getAllOfAnAttribute(stupid_objects,"three")
        self.assertItemsEqual(return_list, range(300))
        
    def test_max_distance(self):
        data = [1, 2, 3, 4, 5,
                   6, 7, 8, 9,
                     10,11,12,
                        20,13,
                           14]
        matrix = CondensedDistanceMatrix(data)
        elements = range(6)
        cluster = Cluster(None,elements)
        self.assertEqual(calculate_max_distance(cluster, matrix),5)

#    def test_generate_separation_report(self):
#        statistics_dic = {}
#        statistics_dic["elems_percent_pure_A"] = 1. 
#        statistics_dic["elems_self_percent_pure_A"] = 2.  
#        statistics_dic["elems_percent_pure_B"] = 3. 
#        statistics_dic["elems_self_percent_pure_B"] = 4. 
#        statistics_dic["elems_percent_mixed"] = 5.
#        s = StatisticsCalculator("matrix","traj_A_pdb","traj_B_pdb", 7, 5)
#        s.statistics_dic= statistics_dic
#        pure1 = [Cluster(None,elements = [3,4,5])]
#        pure2 = [Cluster(None,elements = [8,9])]
#        mixed = [Cluster(None,elements = [0,1,2,7]), Cluster(None,elements = [6,10,11])]
#        expected_string = """
#The clustering has 0 elements (7 from traj_A_pdb and 5 from traj_B_pdb).
#The 1.00% of elements of traj_A_pdb, corresponding to 1 clusters, are separated.
#\t- Which means a 2.00% of its elements.
#The 3.00% of elements of traj_B_pdb, corresponding to 1 clusters, are separated.
#\t- Which means a 4.00% of its elements.
#The 5.00% of elements, corresponding to 2 clusters, are in mixed clusters.
#"""
#        self.assertEqual(expected_string, s.generate_separation_report(pure1,pure2,mixed,0))
#        
#    def test_generate_population_report(self):
#        pass
#    
#    def test_generate_metrics_report(self):
#        pass
    
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_Name']
    unittest.main()