'''
Created on 20/04/2012

@author: victor
'''
import unittest
from pyproct.algorithms.random.RandomAlgorithm import RandomClusteringAlgorithm
from pyproct.matrix.condensedMatrix import CondensedDistanceMatrix
from pyproct.algorithms.random.FakeDistributionRandomAlgorithm import FakeDistributionRandomClusteringAlgorithm


class Test(unittest.TestCase):


    def test_everything(self):
        distances = CondensedDistanceMatrix([ 12.36931688,   5.83095189,   9.43398113,  12.52996409,  15.65247584,
                                             17.4642492,    9.21954446,   4.47213595,   3.16227766,   4.47213595,
                                             5.65685425,   5.,           8.06225775,  11.18033989,  13.15294644,
                                             3.16227766,   6.32455532,   8.24621125,   3.16227766,   5.09901951,   2.  ])
        
        rand_alg = RandomClusteringAlgorithm(distances)
        
        for i in range(100): #@UnusedVariable
            clusterization = rand_alg.perform_clustering(kwargs = {"max_num_of_clusters":5})
            self.assertLess(len(clusterization.clusters),6)
            absolutely_all_elements = []
            for c in clusterization.clusters:
                absolutely_all_elements.extend(c.all_elements)
            absolutely_all_elements.sort()
            self.assertItemsEqual(absolutely_all_elements, range(distances.row_length))
            
    def test_everything_for_rand_distrib(self):
        distances = CondensedDistanceMatrix([0]*int((100*(100-1))/2))
        distribution = [60,33,7]
        rand_alg = FakeDistributionRandomClusteringAlgorithm(distances)
        clusterization  = rand_alg.perform_clustering(kwargs = {"distribution":distribution})
        class_map = {0:0.,1:0.,2:0.}
        
        for i in range(len(clusterization.clusters)):
            class_map[i] += clusterization.clusters[i].get_size()
        
        for i in range(len(clusterization.clusters)):
            self.assertTrue(class_map[i] >distribution[i]-1 and class_map[i]<distribution[i]+1)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_everything']
    unittest.main()