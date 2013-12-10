'''
Created on 23/03/2012

@author: victor
'''
import unittest
from scipy.spatial.distance import pdist
from pyproct.matrix.condensedMatrix import CondensedDistanceMatrix
from pyproct.algorithms.hierarchical.hierarchicalAlgorithm import HierarchicalClusteringAlgorithm
import cStringIO

class HierarchicalStub:
    def __init__(self,cutoff_map):
        self.cutoff_map = cutoff_map

    def perform_clustering(self,kwargs):
        """
        Returns a list with a size determined by the cutoff.
        """
        cutoff = kwargs["cutoff"]
        return ClusteringStub([0]*self.cutoff_map[cutoff])

class ClusteringStub:
    def __init__(self,list):
        self.clusters = list
    def __str__(self):
        return "Clustering ("+str(len(self.clusters))+")"
    def __repr__(self):
        return str(self)

class Test(unittest.TestCase):
    def test_clustering(self):
        """
        Retrocompatibility test.
        MAYBE THIS IMPLEMENTATION OF HCLUSTER DOES NOT WORK THAT WELL
        """
        observations = [(1,1),(1.5,1.5),(5,5),(3,4),(4,4),(3,3.5)]
        matrix_data = pdist(observations)
        condensed = CondensedDistanceMatrix(matrix_data)
        algorithm = HierarchicalClusteringAlgorithm(condensed)
        clusterization = algorithm.perform_clustering(kwargs = {"cutoff":0.5, "method":'single'})
        clusters_string = """[0[0, 1]][3[3, 5]][4[4]][2[2]]"""
        out = cStringIO.StringIO()
        for c in clusterization.clusters:
            out.write( str(c ))
        self.assertEqual(out.getvalue(), clusters_string)



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_clustering']
    unittest.main()