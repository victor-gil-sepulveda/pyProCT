'''
Created on 19/03/2012

@author: victor
'''
import unittest
import pyproclust.algorithms.validation.validationTools as val_tools


class Test(unittest.TestCase):


    def test_dataset_loading_2D(self):
        dataset_string = """41 45
39 44 
42 43
44 43
10 42
"""
        expected_observations = [(41,45),(39,44),(42,43),(44,43),(10,42)]
        observations = val_tools.dataset_loading_2D(dataset_string)
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
        bounding_box = val_tools.get_2D_bounding_box(observations)
        self.assertItemsEqual(expected_bb,bounding_box)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()