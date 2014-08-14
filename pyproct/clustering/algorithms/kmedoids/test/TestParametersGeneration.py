"""
Created on 20/11/2013

@author: victor
"""
import unittest
from pyproct.clustering.algorithms.kmedoids.parametersGeneration import ParametersGenerator


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
        generator  = ParametersGenerator(parameters,{})
        self.assertEqual(2, generator.num_clusters_step)
        self.assertEqual(15, len(generator.get_parameters()[0]))

        del parameters["clustering"]["algorithms"]["kmedoids"]["max"]
        generator  = ParametersGenerator(parameters,{})
        self.assertEqual(2, generator.num_clusters_step)
        self.assertEqual(15, len(generator.get_parameters()[0]))
        self.assertDictEqual({'seeding_type': 'EQUIDISTANT', 'k': 30, 'seeding_max_cutoff': None}, generator.get_parameters()[0][-1])

    def test_random_parameters_generation(self):
        parameters_1 ={
                         "clustering":{
                                       "algorithms":{
                                                     "kmedoids":{
                                                                 "max": 20,
                                                                 "seeding_type":"RANDOM"
                                                                 }
                                                     }
                                       },
                         "evaluation":{
                                       "maximum_clusters": 2,
                                       "minimum_clusters":2
                                       }
                     }

        parameters_2 ={
                         "clustering":{
                                       "algorithms":{
                                                     "kmedoids":{
                                                                 "max": 20,
                                                                 "seeding_type":"RANDOM",
                                                                 "tries":2
                                                                 }
                                                     }
                                       },
                         "evaluation":{
                                       "maximum_clusters":2,
                                       "minimum_clusters":2
                                       }
                     }

        gen_parameters = ParametersGenerator(parameters_1,{}).get_parameters()[0]
        self.assertEqual(10, len(gen_parameters))
        # All seeds are different
        for i in range( 0, len(gen_parameters)-1):
            for j in range( i+1, len(gen_parameters)):
                self.assertNotEqual(gen_parameters[i]["rand_seed"], gen_parameters[j]["rand_seed"])

        gen_parameters = ParametersGenerator(parameters_2,{}).get_parameters()[0]
        self.assertEqual(2, len(gen_parameters))
        # All seeds are different
        for i in range( 0, len(gen_parameters)-1):
            for j in range( i+1, len(gen_parameters)):
                self.assertNotEqual(gen_parameters[i]["rand_seed"], gen_parameters[j]["rand_seed"])



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()