'''
Created on 23/03/2012

@author: victor
'''
import unittest
from scipy.spatial.distance import pdist
from pyproclust.matrix.condensedMatrix import CondensedDistanceMatrix
from pyproclust.algorithms.hierarchical.hierarchicalAlgorithm import HierarchicalClusteringAlgorithm
import cStringIO
from pyproclust.algorithms.hierarchical.hierarchicalTools import calculate_coarse_grain_cutofffs,\
    calculate_fine_grain_cutoffs, dico, calculate_cutoff_numcluster_list
import numpy

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
    
    def test_calculate_coarse_grain_cutofffs(self):
        expected = [ 3.6, 3.2, 2.8, 2.4, 2.0, 1.6, 1.2, 0.8, 0.4, 0]
        numpy.testing.assert_array_almost_equal(expected, calculate_coarse_grain_cutofffs(10, 4),decimal = 3)
    
    def test_calculate_fine_grain_cutoffs(self):
        expected = [2.36, 2.32, 2.28, 2.24, 2.2, 2.16, 2.12, 2.08, 2.04]
        numpy.testing.assert_array_almost_equal(expected, calculate_fine_grain_cutoffs(10,10,4, 2))
        
    def test_dico(self):
        fine_grain_cutoffs = {9:1, 8:1, 7:1, 6:1, 5:1, 4:3, 3:4, 2:10, 1:20, 0:150}
        hierarchical = HierarchicalStub(fine_grain_cutoffs)
        cutoffs = sorted(fine_grain_cutoffs.keys(),reverse = True)
        k = dico(cutoffs, 0, 9, 4, hierarchical, None)
        self.assertEqual(5, k)
        
    def test_calculate_cutoff_numcluster_list(self):
        coarse_grain_cutoffs = {0:200,10:2,20:1,30:1,40:1,50:1,60:1,70:1,80:1,90:1,100:1}
        fine_grain_cutoffs = {19:1, 18:1, 17:1, 16:1, 15:1, 14:3, 13:4, 12:10, 11:20, 10:2}
        all_cutoffs = dict(coarse_grain_cutoffs.items() + fine_grain_cutoffs.items())
        hierarchical = HierarchicalStub(all_cutoffs)
        results = calculate_cutoff_numcluster_list(hierarchical, None, 100, 10, 10)
        expected_cutoffs = [0,10,11,12,13,14]
        cutoffs = []
        for r in results :
            cutoffs .append( r[0])   
        self.assertItemsEqual(expected_cutoffs, cutoffs)
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_clustering']
    unittest.main()