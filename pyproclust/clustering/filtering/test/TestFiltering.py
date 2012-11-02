'''
Created on 07/08/2012

@author: victor
'''
import unittest
from pyproclust.clustering.cluster import Cluster
from pyproclust.clustering.clusterization import Clustering
from pyproclust.protocol.protocolParameters import ProtocolParameters
from pyproclust.clustering.filtering.clusteringFilter import ClusteringFilter

class Test(unittest.TestCase):

    def test_clustering_is_allowed(self):
        clusters =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization1 = Clustering(clusters,details="clustering 1")
        
        clusters =(
                  Cluster(16,[16]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization2 = Clustering(clusters,details="clustering 2")
        
        clusters =(
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization3 = Clustering(clusters,details="clustering 3")
        
        clusters =(
                  Cluster(16,[16,3,5]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization4 = Clustering(clusters,details="clustering 4")
        
        clusters =(
                  Cluster(16,[16,3,5]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization5 = Clustering(clusters,details="clustering 5")
        
        params = ProtocolParameters()
        params.max_clusters = 5
        params.min_clusters = 4
        params.max_noise = 20       
        params.min_cluster_size = 2
        
        cfilter = ClusteringFilter(params)
        
        # non noisy clusters: 3 
        out = cfilter._ClusteringFilter__clusteringIsAllowed(clusterization1, 17)
        self.assertEqual(out, "Out num. clusters range:[4, 5] with 3 clusters\n")
        self.assertEqual(len(clusterization1.clusters), 3)
        
        # non noisy clusters: 2
        expected_out = "Out num. clusters range:[4, 5] with 2 clusters\nToo much noise (max = 20.000%) with 35.294% of noise\n"
        out = cfilter._ClusteringFilter__clusteringIsAllowed(clusterization2, 17)
        self.assertEqual(out, expected_out)
        
        # non noisy clusters: 6, noise  0%
        expected_out ="Out num. clusters range:[4, 5] with 6 clusters\n"
        out = cfilter._ClusteringFilter__clusteringIsAllowed(clusterization3, 31)
        self.assertEqual(out,expected_out)
        
        # non noisy clusters: 4, noise  0%
        out = cfilter._ClusteringFilter__clusteringIsAllowed(clusterization4, 18)
        self.assertEqual(out,"")
        
        # non noisy clusters: 4, noise  19%
        out = cfilter._ClusteringFilter__clusteringIsAllowed(clusterization5, 22)
        self.assertEqual(out,"")
        
        
    def test_clustering_filter(self):
        clusters =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization1 = Clustering(clusters,details="clustering 1")
        
        clusters =(
                  Cluster(16,[16]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization2 = Clustering(clusters,details="clustering 2")
        
        clusters =(
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization3 = Clustering(clusters,details="clustering 3")
        
        clusters =(
                  Cluster(16,[16,3,5]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization4 = Clustering(clusters,details="clustering 4")
        
        clusters =(
                  Cluster(16,[16,3,5]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization5 = Clustering(clusters,details="clustering 5")
        
        params = ProtocolParameters()
        params.max_clusters = 5
        params.min_clusters = 4
        params.max_noise = 20       
        params.min_cluster_size = 2
        
        clusterings = [clusterization1,clusterization2,clusterization3,clusterization4,clusterization5]
        
        cfilter = ClusteringFilter(params)
        filtered = cfilter.doClusteringFiltering(clusterings,17)
        self.assertEqual(len(filtered), 2)
        for c in filtered:
            self.assertIn(c.details, ["clustering 4","clustering 5"])

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()