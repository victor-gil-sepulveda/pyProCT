'''
Created on 17/04/2012

@author: victor
'''
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproclust.algorithms.dbscan.dbscanAlgorithm import DBSCANAlgorithm,\
    PointClassType
# from pyproclust.algorithms.dbscan.dbscanTools import kth_elements_distance, k_dist,\
#     k_scale_gen
from pyproclust.algorithms.dbscan.cython.cythonDbscanTools import kth_elements_distance, k_dist,\
    k_scale_gen
import numpy


class Test(unittest.TestCase):


    def test_eps_neighborhood(self):
        distances = CondensedMatrix([ 1., 1., 1., 1., 1., 1., 1., 1.,
                                         1., 1., 1., 1., 2., 1. ,1.,
                                            1., 1., 1., 3., 1., 1.,
                                               1., 1., 2., 1., 1.,
                                                  1., 1., 1., 1.,
                                                     1., 1., 1.,
                                                        1., 2.,
                                                           1.])
        dbscan_alg = DBSCANAlgorithm(distances)
        
        self.assertItemsEqual(dbscan_alg._DBSCANAlgorithm__eps_neighborhood(6,3),[0, 1, 2, 3, 4, 5, 7, 8])
        self.assertItemsEqual(dbscan_alg._DBSCANAlgorithm__eps_neighborhood(6,2),[0, 1, 3, 4, 5, 7, 8])
        self.assertItemsEqual(dbscan_alg._DBSCANAlgorithm__eps_neighborhood(6,1),[0, 4, 5, 7])
        self.assertItemsEqual(dbscan_alg._DBSCANAlgorithm__eps_neighborhood(6,1.),[0, 4, 5, 7])
        
        self.assertItemsEqual(distances.element_neighbors_within_radius(6,3),[0, 1, 2, 3, 4, 5, 7, 8])
        self.assertItemsEqual(distances.element_neighbors_within_radius(6,2),[0, 1, 3, 4, 5, 7, 8])
        self.assertItemsEqual(distances.element_neighbors_within_radius(6,1),[0, 4, 5, 7])
        self.assertItemsEqual(distances.element_neighbors_within_radius(6,1.),[0, 4, 5, 7])

    def test_seed_expansion(self):
        """
        Graph:
        1--2  3--4
        \ /
         0
         
        Adj. list:
        
        1 1 0 0 0
          1 0 0 0 
            0 0 0 
              0 0 
                1
                
        Inv. distance adjacency list (two nodes are connected if their distance is 
        lower than x)
        
        0 0 x x x
          0 x x x 
            x x x 
              x x 
                0
                      
        """
        distances = CondensedMatrix([ 0., 0., 2., 2., 
                                         0., 2., 2.,
                                            2., 2.,
                                               0.])
        dbscan_alg = DBSCANAlgorithm(distances)
        eps = 1.0
        minpts = 2
        dbscan_alg._DBSCANAlgorithm__seed_expansion(1, eps, minpts, [2])
        expected_classes = [1,1,1,PointClassType.UNCLASSIFIED,PointClassType.UNCLASSIFIED ]
        self.assertItemsEqual(expected_classes, dbscan_alg.element_class)
        
        dbscan_alg.element_class = [PointClassType.UNCLASSIFIED]*5
        dbscan_alg._DBSCANAlgorithm__seed_expansion(2, eps, minpts, [3])
        expected_classes = [PointClassType.UNCLASSIFIED,PointClassType.UNCLASSIFIED,PointClassType.UNCLASSIFIED,PointClassType.UNCLASSIFIED,PointClassType.UNCLASSIFIED ]
        self.assertItemsEqual(expected_classes, dbscan_alg.element_class)
        
        minpts = 1
        dbscan_alg._DBSCANAlgorithm__seed_expansion(2, eps, minpts, [3])
        expected_classes = [PointClassType.UNCLASSIFIED,PointClassType.UNCLASSIFIED,PointClassType.UNCLASSIFIED, 2, 2]
        self.assertItemsEqual(expected_classes, dbscan_alg.element_class)
        
        dbscan_alg.element_class = [PointClassType.UNCLASSIFIED]*5
        dbscan_alg._DBSCANAlgorithm__seed_expansion(1, eps, minpts, [2,4])
        expected_classes = [1, 1, 1, 1, 1]
        self.assertItemsEqual(expected_classes, dbscan_alg.element_class)
        
    def test_dbscan(self):
        distances = CondensedMatrix([ 0., 0., 2., 2., 
                                         0., 2., 2.,
                                            2., 2.,
                                               0.])
        dbscan_alg = DBSCANAlgorithm(distances)
        eps = 1.0
        minpts = 2
        dbscan_alg.perform_clustering(kwargs = {"eps":eps, "minpts":minpts})
        expected = [1, 1, 1, 0, 0]
        self.assertItemsEqual(dbscan_alg.element_class,expected)
        dbscan_alg.element_class = [PointClassType.UNCLASSIFIED]*5
        eps = 1.0
        minpts = 1
        dbscan_alg.perform_clustering(kwargs = {"eps":eps, "minpts":minpts})
        self.assertItemsEqual(dbscan_alg.element_class,[1, 1, 1, 2, 2])
        
    def test_dbscan_regression_mini(self):
        distances = CondensedMatrix([ 12.36931688,   5.83095189,   9.43398113,  12.52996409,  15.65247584,
                                             17.4642492,    9.21954446,   4.47213595,   3.16227766,   4.47213595,
                                             5.65685425,   5.,           8.06225775,  11.18033989,  13.15294644,
                                             3.16227766,   6.32455532,   8.24621125,   3.16227766,   5.09901951,   2.  ])
        dbscan_alg = DBSCANAlgorithm(distances)
        dbscan_alg.perform_clustering(kwargs = {"eps":4.0, "minpts":3})
        self.assertItemsEqual(dbscan_alg.element_class,[0, 1, 0, 1, 1, 1, 0])
        
    def test_kth_elements(self):
        distances = CondensedMatrix([17.46,   9.21,  4.47,  3.16,   4.47,   5.65,   
                                             12.36,  5.83,  9.43,  12.52,  15.65,
                                                     5.,    8.06,  11.18,  13.15,
                                                            3.16,   6.35,   8.24,   
                                                                    3.16,   5.10,
                                                                            2.  ])
        buffer = numpy.empty(distances.row_length)
        expected = [3.1600000000000001, 5.0999999999999996]
        result = kth_elements_distance(4,numpy.array([2,4]),buffer,distances)
        numpy.testing.assert_array_almost_equal(expected, result, decimal = 3)
        
    def test_k_dist(self):
        distances = CondensedMatrix([17.46,   9.21,  4.47,  3.16,   4.47,   5.65,   
                                             12.36,  5.83,  9.43,  12.52,  15.65,
                                                     5.,    8.06,  11.18,  13.15,
                                                            3.16,   6.35,   8.24,   
                                                                    3.16,   5.10,
                                                                            2.  ])
        buffer = numpy.empty(distances.row_length)
        result = k_dist(numpy.array([2,4]), buffer, distances)
        # [[  4.47   5.65]
        # [  9.43  12.52]
        # [  8.06  11.18]
        # [  4.47   5.83]
        # [  3.16   5.1 ]
        # [  3.16   6.35]
        # [  5.1    8.24]]
        expected = [sorted([  4.47,   9.43,  8.06,    4.47,   3.16,   3.16,   5.1]),
                    sorted([  5.65,  12.52, 11.18,    5.83,   5.1,    6.35,   8.24])]
        numpy.testing.assert_array_almost_equal(result,expected,decimal = 2)
        
    def test_k_scale_gen(self):
        self.assertItemsEqual( k_scale_gen(100),    [2, 4, 8, 16, 32, 64])
        self.assertItemsEqual( k_scale_gen(7),      [2, 4])
        self.assertItemsEqual( k_scale_gen(0),      [2])
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testDBSCAN']
    unittest.main()