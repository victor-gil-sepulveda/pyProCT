"""
Created on 13/02/2013

@author: victor
"""
from pyRMSD.matrixHandler import MatrixHandler as pyRMSD_MatrixHandler

class MatrixHandler(object):
    
    def __init__(self, distance_provider, matrix_params):
        """
        Class constructor.
        
        :param distance_provider: A distance provider (something that gives distances by
        indexing a tuple e.g. distance = d_prov[item1,item2].
        
        :param matrix_parameters: The parameters used to build the provider if it is a matrix.
        """
        self.matrix_parameters = matrix_params

        self.distance_matrix = distance_provider

    def save_matrix(self, matrix_path):
        """
        Writes matrix contents to disk.

        @param matrix_save_file: Complete path (with filename) where to save the matrix.
        """
        pyRMSD_MatrixHandler.save_matrix(matrix_path, self.distance_matrix)

    def save_statistics(self, matrix_base_path):
        """
        Writes matrix statistics to disk in JSON format.

        @param matrix_base_path: The folder where to save the 'statistics.json' file.
        """
        return pyRMSD_MatrixHandler.save_statistics(matrix_base_path, self.distance_matrix)
    

