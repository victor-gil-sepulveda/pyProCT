"""
Created on 08/01/2013

@author: victor
"""
import time
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.clustering.cluster import Cluster
from pyproct.clustering.clustering import Clustering
from  pyproct.clustering.metrics.cython.normNCut import CythonNCut
from pyproct.clustering.metrics.cython.boundedCohesion import CythonBoundedCohesionCalculator
from pyproct.clustering.metrics.cython.meanMinimumDistance import CythonMeanMinimumDistanceCalculator
from pyproct.clustering.metrics.cython.silhouette import CythonSilhouetteCoefficientCalculator
from pyproct.clustering.metrics.graphMetrics import NCut
from pyproct.clustering.metrics.boundedCohesion import BoundedCohesionCalculator
from pyproct.clustering.metrics.meanMinimumDistance import MeanMinimumDistanceCalculator
from pyproct.clustering.metrics.silhouette import SilhouetteCoefficientCalculator


def do_benchmark( who, calculator, matrix, clustering):
    t1 = time.time()
    result = calculator.evaluate(clustering, matrix)
    t2 = time.time()
    print '\tCalculating with %s took %0.3fs giving %.3f as result' % (who,t2-t1,result)
    return t2-t1

if __name__ == '__main__':
    
    # Create a random matrix
    MAX_ELEMENTS = 10000
    DATA_LEN = (MAX_ELEMENTS * (MAX_ELEMENTS-1))/2
    matrix = CondensedMatrix(numpy.random.rand(DATA_LEN))
    
    #Create a random clustering
    NUMBER_OF_CLUSTERS = 100
    ELEMENTS_PER_CLUSTER = MAX_ELEMENTS/ NUMBER_OF_CLUSTERS
    clusters = []
    for i in range(NUMBER_OF_CLUSTERS):
        elems = range(i*ELEMENTS_PER_CLUSTER,(i+1)*ELEMENTS_PER_CLUSTER)
        clusters.append(Cluster(None, elems))
    clustering = Clustering(clusters, "Fake clustering")
    
    # Create metrics
    metrics = {
               "NCut":(CythonNCut(),NCut()),
               "BoundedCohesion":(CythonBoundedCohesionCalculator(),BoundedCohesionCalculator()),
               "MeanMinimumDistance":(CythonMeanMinimumDistanceCalculator(), MeanMinimumDistanceCalculator()),
               "Silhouette":(CythonSilhouetteCoefficientCalculator(),SilhouetteCoefficientCalculator())
               }
    
    # Do the benchmark
    for key in metrics:
        print "Benchmarking "+key+":"
        cythoncalc, calc = metrics[key]
        cytime = do_benchmark("Cython::"+key, cythoncalc, matrix, clustering)
        pytime = do_benchmark("Python::"+key,calc, matrix, clustering)
        print "\tSpeedup: %.1fX"%(pytime/cytime)
        