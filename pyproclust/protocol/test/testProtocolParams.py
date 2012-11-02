'''
Created on 04/06/2012

@author: victor
'''
import unittest
from pyproclust.protocol.protocolParameters import ProtocolParameters
from pyproclust.clustering.clusterization import Clustering
from pyproclust.clustering.cluster import Cluster


class Test(unittest.TestCase):


    def test_choose_best_clustering(self):
        c1 = "c1"
        c2 = "c2"
        c3 = "c3"
        results_pack = [
                        (c1,{"metric1":0.5, "metric2":0.2, "metric3":0.2}),
                        (c2,{"metric1":0.33, "metric2": 0.33, "metric3": 0.33}),
                        (c3,{"metric1":0.2, "metric2": 0.3, "metric3": 0.5}),
                        ]
        
        params = ProtocolParameters()
        params.cluster_score_value_map = {"metric1":1, "metric2":1, "metric3":1}
        self.assertEqual( "c3",params.chooseBestClustering(results_pack)[1])
        params.cluster_score_value_map = {"metric1":0.5, "metric2":1, "metric3":0}
        self.assertEqual( "c2",params.chooseBestClustering(results_pack)[1])
        params.cluster_score_value_map = {"metric1":0, "metric2":1, "metric3":1}
        self.assertEqual( "c3",params.chooseBestClustering(results_pack)[1])
        
        
         
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_chooseBestClustering']
    unittest.main()