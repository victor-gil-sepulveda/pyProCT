'''
Created on 19/11/2013

@author: victor
'''
import os
import prody
import unittest
import numpy
from pyproclust.clustering.cluster import Cluster
from pyproclust.tools.test.data.amber_short_ca import amber_short_ca_contents
from pyproclust.clustering.comparison.caDisplacement import CA_mean_square_displacement_of_cluster,\
    calc_norms

class Test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # Generate and read the pdb
        cls.pdb_path = "tmp_pdb.pdb"
        open(cls.pdb_path,"w").write(amber_short_ca_contents);
        cls.pdb = prody.parsePDB(cls.pdb_path)
        os.system("rm %s"%cls.pdb_path)

    def test_CA_mean_square_displacement_of_cluster(self):
        alpha_carbons_trajectory_pdb = self.pdb.select("name CA")
        expected = [ 1.32850684, 1.56350125, 1.51248423, 1.96403403, 1.81467535, 2.44374607,
                     3.36425126, 3.95337708, 3.3888146,  3.10514178, 2.69929018, 2.35866313,
                     2.23312347, 2.05268284, 1.55179237, 1.3045629,  0.66721661, 0.33162028,
                     0.8483944,  1.26908385, 1.16434935, 1.69391873, 1.36914871, 2.12399383,
                     2.15054364, 1.54744313, 2.00120183, 2.68321728, 2.27389936, 2.59041278,
                     3.46147909, 3.71143747, 3.59721538, 3.80392142, 4.1446817,  3.61275165,
                     3.68878658, 3.01804669, 2.95709462, 2.68221771, 2.28039218, 1.61026319,
                     1.16795743, 1.70812982, 2.03046519, 3.17433966, 3.35453324, 2.69025038,
                     1.82128602, 1.39444184, 1.66178521, 1.86392499, 2.51830867, 2.04003363,
                     1.6348583,  1.37095361, 1.84063533, 2.02663961, 1.44658883, 1.87954519,
                     1.51140317, 1.7199428,  1.71083824, 2.161764,   2.305442,   2.34243684,
                     1.64746632, 2.01549208, 1.8223211,  2.20481542, 2.53409383, 2.38280986,
                     3.19762438, 3.08093492, 3.17642258, 4.22996496]
        numpy.testing.assert_array_almost_equal(expected, CA_mean_square_displacement_of_cluster(alpha_carbons_trajectory_pdb, Cluster(2,[2,3,5])))
    
    def test_calc_norms(self):
        cooordsets = numpy.array([[[1,2,3],
                                   [4,5,6],
                                   [7,8,9]],
                                  [[10,11,12],
                                   [13,14,15],
                                   [16,17,18]],
                                  [[19,20,21],
                                   [22,23,24],
                                   [25,26,27]],
                                  [[19,20,21],
                                   [22,23,24],
                                   [25,26,27]]])
        expected = [ 19.48557159, 19.48557159, 19.48557159]
        numpy.testing.assert_array_almost_equal(expected, calc_norms(cooordsets,cooordsets[0]))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()