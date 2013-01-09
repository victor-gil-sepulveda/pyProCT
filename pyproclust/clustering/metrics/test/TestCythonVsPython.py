'''
Created on 09/01/2013

@author: victor
'''
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
from  pyproclust.clustering.metrics.cython.normNCut import CythonNCut
from pyproclust.clustering.metrics.cython.boundedCohesion import CythonBoundedCohesionCalculator
from pyproclust.clustering.metrics.cython.meanMinimumDistance import CythonMeanMinimumDistanceCalculator
from pyproclust.clustering.metrics.cython.silhouette import CythonSilhouetteCoefficientCalculator
from pyproclust.clustering.metrics.graphMetrics import NCut
from pyproclust.clustering.metrics.boundedCohesion import BoundedCohesionCalculator
from pyproclust.clustering.metrics.meanMinimumDistance import MeanMinimumDistanceCalculator
from pyproclust.clustering.metrics.silhouette import SilhouetteCoefficientCalculator
from pyproclust.clustering.clustering import Clustering
from pyproclust.clustering.cluster import Cluster
import numpy

class TestCythonVsPython(unittest.TestCase):
    """ These tests suppose that the python code is tested and working """

    def test_cython_vs_python(self):
        #metrics with a Cython implementation
        metrics = {
               "NCut":(CythonNCut(),NCut()),
               "BoundedCohesion":(CythonBoundedCohesionCalculator(),BoundedCohesionCalculator()),
               #"MeanMinimumDistance":(CythonMeanMinimumDistanceCalculator(), MeanMinimumDistanceCalculator()),
               "Silhouette":(CythonSilhouetteCoefficientCalculator(),SilhouetteCoefficientCalculator())}
        
        for i in range(20):
            print "\n",i, "try"
            matrix = CondensedMatrix(numpy.random.rand(1000))
            clusters = []
            for i in range(10):
                elems = range(i*100,(i+1)*100)
                clusters.append(Cluster(None, elems))
            clustering = Clustering(clusters, "Fake clustering")
            for metric_name in metrics:
                print "\t"+metric_name
                cythoncalc, calc = metrics[metric_name]
                cyresult = cythoncalc.evaluate(clustering, matrix)
                result = calc.evaluate(clustering, matrix)
                self.assertEqual(cyresult, result)
                

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()