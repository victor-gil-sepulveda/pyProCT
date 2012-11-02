'''
Created on 07/09/2012

@author: victor
'''
import unittest
from pyproclust.clustering.selection.bestClusteringSelector import BestClusteringSelector


class Test(unittest.TestCase):

    def test_init(self):
        value_map = {"A":(2.,">"), "B":(1.0,"<")}
        sel = BestClusteringSelector([value_map])
        self.assertEqual(sel.calcMaxScoring(value_map),3.)
    
    def test_score_clustering(self):
        value_map  = {"A":(2.,">"), "B":(1.0,"<")}
        sel = BestClusteringSelector([value_map])
        # Eval list for clustering
        self.assertEqual( sel.scoreClustering({"A":0.5,"B":0.7,"C":0.3},value_map), (0.5*2 + (1-0.7)*1)/3)
    
    def test_get_best_clustering(self):
        value_map  = {"A":(2.,">"), "B":(1.0,"<")}
        sel = BestClusteringSelector([value_map])
        results_pack = [("Clustering 1", {"A":0.5,"B":0.7,"C":0.3}),
                        ("Clustering 2", {"A":0.7,"B":0.3,"C":0.3}),
                        ("Clustering 3", {"A":0.7,"B":0.7,"C":0.3})]
        best_score, best_clustering =  sel.chooseBestClustering(results_pack) #@UnusedVariable
        self.assertEqual("Clustering 2",best_clustering)
        
        value_maps  = [{"A":(2.,">"), "B":(1.0,"<")},{"C":(2.0,">")}]
        sel = BestClusteringSelector(value_maps)
        results_pack = [("Clustering 1", {"A":0.5,"B":0.7,"C":0.3}),
                        ("Clustering 2", {"A":0.7,"B":0.3,"C":0.3}),
                        ("Clustering 3", {"A":0.7,"B":0.7,"C":1.0})]
        best_score, best_clustering =  sel.chooseBestClustering(results_pack) #@UnusedVariable
        self.assertEqual("Clustering 3",best_clustering)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()