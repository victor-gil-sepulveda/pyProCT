"""
Created on 27/01/2014

@author: victor
"""
import unittest

from pyproct.clustering.clustering import Clustering
import pyproct.clustering.metrics.test.data.example_clustering_1 as data
from pyproct.driver.handlers.matrix.matrixHandler import MatrixHandler
from pyproct.clustering.metrics.silhouette import SilhouetteCoefficientCalculator


class Test(unittest.TestCase):

    def testSilhouetteSpecialCase(self):
        clustering = Clustering.from_dic(data.clustering_01)
        mh = MatrixHandler({
                                "method": "load",
                                "parameters":{
                                    "path": "data/example_clustering_1_matrix"
                                }
                            }
        )
        s = SilhouetteCoefficientCalculator()
        matrix =  mh.create_matrix(None)
        print s.evaluate(clustering, matrix)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSilhouetteSpecialCase']
    unittest.main()