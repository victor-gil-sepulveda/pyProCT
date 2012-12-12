'''
Created on 31/01/2012

@author: victor
'''
import unittest

import imageStub
from pyproclust.matrix.completeMatrix import CompleteDistanceMatrix
from pyproclust.matrix.condensedMatrix import CondensedDistanceMatrix
from pyproclust.clustering.plotting.matrixPlot import condensed_matrix_to_image,\
    data_matrix_to_image, hvs_to_rgb_converter_yellow_to_red

class Test(unittest.TestCase):

    def __test_color_descriptor(self,norm_value):
        return norm_value
        
    def test_condensed_matrix_to_image(self):
        expected_image_matrix = CompleteDistanceMatrix([[0.0, 1.0, 4.5, 7.2, 6.7], 
                                [1.0, 0.0, 8.5, 4.5, 3.6], 
                                [4.5, 8.5, 0.0, 7.8, 2.2], 
                                [7.2, 4.5, 7.8, 0.0, 2.0], 
                                [6.7, 3.6, 2.2, 2.0, 0.0]])
        
        condensed_matrix = CondensedDistanceMatrix([1.0, 4.5, 7.2, 6.7, 
                                 8.5, 4.5, 3.6, 
                                      7.8, 2.2, 
                                           2.0],True) 
        
        image_pixels = imageStub.ImageStub(5,5)
        
        condensed_matrix_to_image(condensed_matrix, image_pixels, self.__test_color_descriptor)
        
        self.assertEqual(expected_image_matrix,CompleteDistanceMatrix(image_pixels.matrix))

    def test_matrix_to_image(self):
        expected_image_matrix = CompleteDistanceMatrix([[0.0, 1.0, 4.5, 7.2, 6.7], 
                                [1.0, 0.0, 8.5, 4.5, 3.6], 
                                [4.5, 8.5, 0.0, 7.8, 2.2], 
                                [7.2, 4.5, 7.8, 0.0, 2.0], 
                                [6.7, 3.6, 2.2, 2.0, 0.0]])
        
        source_matrix = [[0.0, 1.0, 4.5, 7.2, 6.7], 
                         [1.0, 0.0, 8.5, 4.5, 3.6], 
                         [4.5, 8.5, 0.0, 7.8, 2.2], 
                         [7.2, 4.5, 7.8, 0.0, 2.0], 
                         [6.7, 3.6, 2.2, 2.0, 0.0]]
        
        image_pixels = imageStub.ImageStub(5,5)
        
        data_matrix_to_image(source_matrix, image_pixels, self.__test_color_descriptor)
    
        self.assertEqual(expected_image_matrix, CompleteDistanceMatrix(image_pixels.matrix))

    def test_hvs_to_rgb_converter_yellow_to_red(self):
        one = hvs_to_rgb_converter_yellow_to_red(0.2)
        self.assertEqual( one, (204, 65, 40) )
        two = hvs_to_rgb_converter_yellow_to_red(0.4)
        self.assertEqual( two, (204, 89, 40) )
        three = hvs_to_rgb_converter_yellow_to_red(0.7)
        self.assertEqual( three, (204, 126, 40) )

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()