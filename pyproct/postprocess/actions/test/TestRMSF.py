"""
Created on 19/11/2013

@author: victor
"""
import os
import prody
import unittest
import numpy
from pyproct.clustering.cluster import Cluster
from pyproct.tools.test.data.amber_short_ca import amber_short_ca_contents,\
    amber_short_ca_2
import cStringIO
from pyproct.postprocess.actions.rmsf import superpose_and_calc_rmsf

class TestRMSF(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Generate and read the pdb
        cls.pdb_path = "tmp_pdb.pdb"
        open(cls.pdb_path,"w").write(amber_short_ca_contents);
        cls.pdb = prody.parsePDB(cls.pdb_path)
        os.system("rm %s"%cls.pdb_path)

    def test_rmsf_of_cluster(self):
        """
        Regression test. Can be checked against Prody. 
        """
#         stream = cStringIO.StringIO(amber_short_ca_2)
#         pdb =prody.parsePDBStream(stream)
#         ensemble = prody.Ensemble("")
#         ensemble.setCoords(pdb.getCoords() )
#         ensemble.addCoordset(pdb.getCoordsets())
#         ensemble.iterpose()
#         expected = prody.measure.calcRMSF(ensemble)
        
        expected = [0.641553048, 0.44227368, 0.21570513, 0.18929423, 
                    0.240815418, 0.40685254, 0.49579797, 0.91562725,
                    0.807144101, 1.06956639, 0.72968014, 0.58836750,
                    0.575403705, 0.32143319, 0.28998106, 0.35767000,
                    0.239870647, 0.44757687, 0.50674987, 0.44649181,
                    0.441386524, 0.49127174, 0.19687866, 0.28049732,
                    0.274612892, 0.24545861, 0.26809383, 0.32588850,
                    0.315689301, 0.45270097, 0.62447889, 0.76968487,
                    0.649044680, 0.41381866, 0.54485095, 0.40434944,
                    0.496704737, 0.45709734, 0.39371106, 0.43722143,
                    0.423775698, 0.43673013, 0.32709396, 0.33004627,
                    0.395588502, 0.64509831, 0.92129826, 0.76995881,
                    0.332302544, 0.30982020, 0.30092531, 0.20037805,
                    0.491437813, 0.29777384, 0.26834788, 0.33844817,
                    0.357405734, 0.63051909, 0.60519168, 0.33503306,
                    0.300343918, 0.30988569, 0.22061281, 0.31602567,
                    0.225566485, 0.29537637, 0.28711530, 0.30667463,
                    0.473414725, 0.46727466, 0.49955874, 0.58760692,
                    1.069410813, 1.42671321, 1.92347804, 2.14886810 ]

        
        alpha_carbons_trajectory_pdb = self.pdb.select("name CA").getCoordsets()
        
        numpy.testing.assert_array_almost_equal(expected,
                                                superpose_and_calc_rmsf(alpha_carbons_trajectory_pdb, 
                                                                        Cluster(2,[2,3,5])),
                                                8)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
