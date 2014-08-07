"""
Created on 23/03/2012

@author: victor
"""
import unittest
from scipy.spatial.distance import pdist
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.algorithms.hierarchical.hierarchicalAlgorithm import HierarchicalClusteringAlgorithm
import cStringIO
from pyproct.algorithms.hierarchical.hierarchicalTools import get_cutoff_range
import math

class HierarchicalStub:
    def __init__(self,cutoff_map):
        self.cutoff_map = cutoff_map

    def perform_clustering(self, kwargs):
        """
        Returns a list with a size determined by the cutoff.
        """
        cutoff = kwargs["cutoff"]
        return ClusteringStub([0]*self.cutoff_map[math.ceil(cutoff)])

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
        condensed = CondensedMatrix(matrix_data)
        algorithm = HierarchicalClusteringAlgorithm(condensed)
        clusterization = algorithm.perform_clustering(kwargs = {"cutoff":0.5, "method":'single'})
        clusters_string = """[0[0, 1]][3[3, 5]][4[4]][2[2]]"""
        out = cStringIO.StringIO()
        for c in clusterization.clusters:
            out.write(str(c))
        self.assertEqual(out.getvalue(), clusters_string)

    def test_cutoff_range(self):
        fine_grain_cutoffs = {9:1, 8:1, 7:1, 6:1, 5:1, 4:3, 3:4, 2:10, 1:20, 0:150}
        hierarchical = HierarchicalStub(fine_grain_cutoffs)
        (left,right) = get_cutoff_range(starting_cutoff = 0,
                             ending_cutoff = 9,
                             min_clusters = 3,
                             max_clusters = 7,
                             grain = 0.01,
                             hie_algorithm=hierarchical)
        expected = (2,4)
        self.assertEqual(expected,(math.floor(left),math.ceil(right)))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_clustering']
    unittest.main()