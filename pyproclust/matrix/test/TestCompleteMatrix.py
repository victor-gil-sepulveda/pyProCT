'''
Created on 27/01/2012

@author: victor
'''
import unittest
import cStringIO
from pyproclust.matrix.completeMatrix import CompleteDistanceMatrix,load_complete_matrix
from pyproclust.matrix.matrixCommon import assertMatrixAlmostEqual


class TestMatrixOperations(unittest.TestCase):

    def test_equal(self):
        matrix1 = CompleteDistanceMatrix([[1,2,3,4,5],
                                          [6,4,0,4,2],
                                          [7,5,1,4,7],
                                          [8,7,5,0,9]])
        matrix2 = CompleteDistanceMatrix([[1,2,3,4,5],
                                          [6,4,0,4,2],
                                          [7,5,1,4,7],
                                          [8,7,5,0,9]])
        reflected = CompleteDistanceMatrix([[8,7,5,0,9],
                                            [7,5,1,4,7],
                                            [6,4,0,4,2],
                                            [1,2,3,4,5]])
        self.assertEqual(matrix1 == matrix2,  True)
        self.assertEqual(matrix2 == matrix1,  True)
        self.assertEqual(matrix1 != reflected,True)
        
    def test_reflect(self):
        matrix = CompleteDistanceMatrix([[1,2,3,4,5],
                                          [6,4,0,4,2],
                                          [7,5,1,4,7],
                                          [8,7,5,0,9]])
        expected = CompleteDistanceMatrix([[8,7,5,0,9],
                                            [7,5,1,4,7],
                                            [6,4,0,4,2],
                                            [1,2,3,4,5]])
        
        matrix.reflect()
        self.assertEqual(matrix == expected,True)
        
    def test_validate_distance_matrix(self):
        matrix1 = CompleteDistanceMatrix([[0.0, 1.0, 5.0, 8.5, 7.2], 
                    [1.0, 0.0, 4.5, 7.8, 6.7], 
                    [5.0, 4.5, 0.0, 3.6, 2.2], 
                    [8.5, 7.8, 3.6, 0.0, 2.0], 
                    [7.2, 6.7, 2.2, 2.0, 0.0]])
        
        matrix2 = CompleteDistanceMatrix([[0, 1, 2, 3, 4],
                   [1, 0, 2, 3, 4],
                   [1, 0, 0, 3, 4],
                   [1, 4, 2, 0, 4],
                   [1, 6, 2, 3, 0]])
        
        self.assertEqual(True,matrix1.symmetry_error()< 0.01)
        self.assertEqual(True,matrix2.symmetry_error() > 1.5)
  
    def test_validate_dimensions(self):
        
        matrix = CompleteDistanceMatrix([[1,2,3,4,5],
                                          [6,4,0,4,2],
                                          [7,5,1,4,7],
                                          [8,7,5,0,9]])
        dim = (len(matrix.get_data()[0]),matrix.get_row_dimension())
        self.assertEqual(True,matrix.validate_dimensions(dim[0],dim[1],verbose = False))
        self.assertEqual(False,matrix.validate_dimensions(3,6,verbose = False))
        
    def test_minmax(self):
        matrix = CompleteDistanceMatrix([[1,2,3,4,5],
                                          [6,4,5,4,2],
                                          [7,5,1,4,7],
                                          [8,7,5,6,9]])
        
        expected = (1,9)
        
        self.assertEqual(matrix.get_minimum_and_maximum(),expected)
    
    def test_normalize_matrix(self):
        matrix = CompleteDistanceMatrix([[1,2,3,4,5],
                                      [6,7,8,9,10],
                                      [11,12,13,14,15],
                                      [16,17,18,19,20]])
        
        expected_matrix = CompleteDistanceMatrix([[0.0, 0.052631578947368418, 0.10526315789473684, 0.15789473684210525, 0.21052631578947367], 
                           [0.26315789473684209, 0.31578947368421051, 0.36842105263157893, 0.42105263157894735, 0.47368421052631576], 
                           [0.52631578947368418, 0.57894736842105265, 0.63157894736842102, 0.68421052631578949, 0.73684210526315785], 
                           [0.78947368421052633, 0.84210526315789469, 0.89473684210526316, 0.94736842105263153, 1.0]])
        
        matrix.normalize((1,20))
        
        assertMatrixAlmostEqual(self,matrix.get_data(), expected_matrix.get_data(), 4)
    
    def test_save_matrix(self):
        matrix_to_save = CompleteDistanceMatrix([[0.0, 1.0, 5.0, 8.5, 7.2], 
                                                [1.0, 0.0, 4.5, 7.8, 6.7], 
                                                [5.0, 4.5, 0.0, 3.6, 2.2], 
                                                [8.5, 7.8, 3.6, 0.0, 2.0], 
                                                [7.2, 6.7, 2.2, 2.0, 0.0]])
        
        # with final spaces!
        expected_string =  """0.0 1.0 5.0 8.5 7.2 
1.0 0.0 4.5 7.8 6.7 
5.0 4.5 0.0 3.6 2.2 
8.5 7.8 3.6 0.0 2.0 
7.2 6.7 2.2 2.0 0.0 
"""
        output = cStringIO.StringIO()
        
        matrix_to_save.save(output)
        
        self.assertEqual(expected_string,output.getvalue())
        
    def test_load_complete_matrix(self):
        matrix_string =  """0.0 1.0 5.0 8.5 7.2
1.0 0.0 4.5 7.8 6.7
5.0 4.5 0.0 3.6 2.2
8.5 7.8 3.6 0.0 2.0
7.2 6.7 2.2 2.0 0.0
"""
        expected = CompleteDistanceMatrix([[0.0, 1.0, 5.0, 8.5, 7.2], 
                                            [1.0, 0.0, 4.5, 7.8, 6.7], 
                                            [5.0, 4.5, 0.0, 3.6, 2.2], 
                                            [8.5, 7.8, 3.6, 0.0, 2.0], 
                                            [7.2, 6.7, 2.2, 2.0, 0.0]])
        input = cStringIO.StringIO(matrix_string)
        matrix = load_complete_matrix(input)
        assertMatrixAlmostEqual(self,matrix.get_data(), expected.get_data(),1)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()