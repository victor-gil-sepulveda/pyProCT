"""
Created on 17/09/2012

@author: victor
"""
import unittest
from pyproct.protocol.refinement.Refiner import Refiner
from pyproct.clustering.cluster import Cluster
from pyproct.clustering.clustering import Clustering

class ClusterStub():
    def __init__(self,elements):
        self.all_elements = elements

class KMedoidsAlgorithmStub():
    def __init__(self, submatrix):
        pass

    def perform_clustering(self, params):
        return Clustering([Cluster(None,[0,1,4]),Cluster(None,[2,3])])

class Test(unittest.TestCase):

    def test_redefine_cluster_with_map(self):
        initial_cluster = Cluster(None,[1,3,4,7,8])
        final_cluster_1 = Cluster(None,[0,1,4]) #-> elements [1,3,8] of initial cluster
        final_cluster_2 = Cluster(None,[2,3]) #-> elements [4,7] of initial cluster

        self.assertItemsEqual( [1,3,8],Refiner.redefine_cluster_with_map(initial_cluster, final_cluster_1).all_elements)
        self.assertItemsEqual( [4,7],Refiner.redefine_cluster_with_map(initial_cluster, final_cluster_2).all_elements)

    def test_repartition_with_kmedoids(self):
        Refiner.KMedoidsAlgorithmClass = KMedoidsAlgorithmStub
        clustering = Refiner.repartition_with_kmedoids(Cluster(None,[1,3,4,7,8]), 0, None)
        # We'll suppose that the order is always the same
        self.assertItemsEqual( [1,3,8],clustering.clusters[0])
        self.assertItemsEqual( [4,7],clustering.clusters[1])



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_Name']
    unittest.main()