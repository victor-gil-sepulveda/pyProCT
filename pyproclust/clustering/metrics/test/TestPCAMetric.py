'''
Created on 20/12/2012

@author: victor
'''
import unittest
import numpy
import prody

from pyRMSD.utils.proteinReading import flattenCoords

from pyproclust.tools.test.data.amber_short_ca import amber_short_ca_contents
from pyproclust.clustering.metrics.pcaMetrics import PCAMetric
from pyproclust.clustering.cluster import Cluster
from pyproclust.clustering.clustering import Clustering
import os

class TrajectoryHandlerStub:
    def __init__(self, coords, apc):
        self.coordinates = coords
        self.atoms_per_conformation = apc

class testPCAMetric(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Generate and read the pdb
        cls.pdb_path = "tmp_pdb.pdb"
        open(cls.pdb_path,"w").write(amber_short_ca_contents);
        try:
            prody.setVerbosity('none')
        except Exception :
            pass 
        cls.pdb = prody.parsePDB(cls.pdb_path, subset='calpha')
        
        # Save coordsets before superposition
        cls.not_iterposed_coordsets = numpy.array(cls.pdb.getCoordsets())
        
        # Do Prody iterposition
        cls.ensemble = prody.Ensemble('pca_test_ensemble')
        cls.ensemble.setCoords( cls.pdb.getCoords())
        cls.ensemble.addCoordset(cls.pdb.getCoordsets())
        #prody.setVerbosity('info')
        cls.ensemble.iterpose()
        cls.coordsets = cls.ensemble.getCoordsets()
        
    @classmethod
    def tearDownClass(cls):
        os.system("rm "+cls.pdb_path)
        
    def test_covariance_matrix_vs_prody(self):
        # do it with PCA metric
        my_cov_matrix = PCAMetric.create_covariance_matrix(testPCAMetric.coordsets)
         
        # Do it with prody
        pca = prody.PCA('pcametric_pca')
        pca.buildCovariance(testPCAMetric.ensemble)
        prody_cov_matrix = pca._cov
         
        # Compare
        numpy.testing.assert_almost_equal(my_cov_matrix, prody_cov_matrix,10)
     
    def test_eigenvalues(self):
        # do it with PCA metric
        my_cov_matrix = PCAMetric.create_covariance_matrix(testPCAMetric.coordsets)
        biggest_eigenvalue = PCAMetric.calculate_biggest_eigenvalue(my_cov_matrix)
         
        # Do it with prody
        pca = prody.PCA('pcametric_pca')
        pca.buildCovariance(testPCAMetric.ensemble)
        pca.calcModes(n_modes=1)
        self.assertAlmostEqual(pca.getEigvals()[0], biggest_eigenvalue,10)
         
    def test_iterative_superposition(self):
        fcoords = flattenCoords(testPCAMetric.not_iterposed_coordsets)
        PCAMetric.do_iterative_superposition(fcoords, 12, 66)
        res_fcoords = fcoords.reshape((12,66,3))
        numpy.testing.assert_almost_equal(res_fcoords, testPCAMetric.coordsets,12)
         
    def test_PCA(self):
        """
        Regression test.
        """
        trajectory_handler = TrajectoryHandlerStub(testPCAMetric.not_iterposed_coordsets,66)
        clustering = Clustering([Cluster(None,range(6)),Cluster(None,range(6,12))], "a clustering")
        pcaMetric = PCAMetric()
        self.assertAlmostEquals(pcaMetric.evaluate(clustering, trajectory_handler), 1.42774868808, 10)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()