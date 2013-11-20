'''
Created on 20/11/2013

@author: victor
'''
import unittest
from pyproclust.algorithms.kmedoids.parametersGeneration import ParametersGenerator


class Test(unittest.TestCase):


    def test_parameters_generation(self):
        parameters ={
                     "clustering":{
                                   "algorithms":{
                                                 "kmedoids":{
                                                             "max": 20 
                                                             }
                                                 }
                                   },
                     "evaluation":{
                                   "maximum_clusters": 30,
                                   "minimum_clusters":2
                                   }
                     }
        generator  = ParametersGenerator(parameters,{});
        self.assertEqual(2, generator.num_clusters_step)
        self.assertEqual(15, len(generator.get_parameters()[0]))
        
        del parameters["clustering"]["algorithms"]["kmedoids"]["max"]
        generator  = ParametersGenerator(parameters,{});
        self.assertEqual(2, generator.num_clusters_step)
        self.assertEqual(15, len(generator.get_parameters()[0]))
        self.assertDictEqual({'seeding_type': 'EQUIDISTANT', 'k': 30, 'seeding_max_cutoff': None}, generator.get_parameters()[0][-1])


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()