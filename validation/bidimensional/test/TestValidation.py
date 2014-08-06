"""
Created on 19/03/2012

@author: victor
"""
import unittest
from validation.validationTools import params_to_string, dataset_loading_2D,\
    get_2D_bounding_box


class Test(unittest.TestCase):


    def test_dataset_loading_2D(self):
        dataset_string = """41 45
39 44 
42 43
44 43
10 42
"""
        expected_observations = [(41,45),(39,44),(42,43),(44,43),(10,42)]
        observations = dataset_loading_2D(dataset_string)
        self.assertItemsEqual(expected_observations, observations)
        
    def test_get_2D_bounding_box(self):
        observations = (
                        (3,4),
                        (5,7),
                        (2,1),
                        (6,0),
                        (43,3),
                        (3,6),
                        (5,10)
                        )
        expected_bb = ((2, 0), (43, 10))
        bounding_box = get_2D_bounding_box(observations)
        self.assertItemsEqual(expected_bb,bounding_box)
        
    def test_params_to_string(self):
        params = {'k':3,"lol":10.,'type':'mytype'}
        self.assertEqual("k_3_type_mytype_lol_10.000",  params_to_string(params))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()