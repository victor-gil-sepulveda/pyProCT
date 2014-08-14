"""
Created on 09/01/2013

@author: victor
"""
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
from  pyproct.clustering.evaluation.metrics.cython.normNCut import CythonNCut
from pyproct.clustering.evaluation.metrics.cython.boundedCohesion import CythonMirrorCohesionCalculator
from pyproct.clustering.evaluation.metrics.cython.meanMinimumDistance import CythonMeanMinimumDistanceCalculator
from pyproct.clustering.evaluation.metrics.cython.silhouette import CythonSilhouetteCoefficientCalculator
from pyproct.clustering.evaluation.metrics.graphMetrics import NCut
from pyproct.clustering.evaluation.metrics.boundedCohesion import MirrorCohesionCalculator
from pyproct.clustering.evaluation.metrics.meanMinimumDistance import MeanMinimumDistanceCalculator
from pyproct.clustering.evaluation.metrics.silhouette import SilhouetteCoefficientCalculator
from pyproct.clustering.clustering import Clustering
from pyproct.clustering.cluster import Cluster
import numpy

class TestCythonVsPython(unittest.TestCase):
    """ These tests suppose that the python code is tested and working """

    def test_cython_vs_python(self):
        #metrics with a Cython implementation
        metrics = {
               "NCut":(CythonNCut(),NCut()),
               "BoundedCohesion":(CythonMirrorCohesionCalculator(),MirrorCohesionCalculator()),
               "MeanMinimumDistance":(CythonMeanMinimumDistanceCalculator(10), MeanMinimumDistanceCalculator(10)),
               "Silhouette":(CythonSilhouetteCoefficientCalculator(),SilhouetteCoefficientCalculator())}
        
        for i in range(5):
            print "\n%dth try"%i
            matrix = CondensedMatrix(numpy.random.rand((1000*999)/2))
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
                if metric_name == "MeanMinimumDistance":
                    # Use of randomness can change results a bit
                    self.assertAlmostEqual(cyresult, result, 2)
                else:
                    self.assertAlmostEqual(cyresult, result, 8)
                

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()