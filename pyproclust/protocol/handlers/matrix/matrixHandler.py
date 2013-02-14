'''
Created on 13/02/2013

@author: victor
'''
import pyRMSD.matrixHandler.MatrixHandler as pyRMSD_MatrixHandler
from pyproclust.protocol.handlers.matrix.rmsdMatrixBuilder import RMSDMatrixBuilder
from pyproclust.protocol.handlers.matrix.euclideanDistanceMatrixBuilder import EuclideanDistanceMatrixBuilder

class MatrixHandler(object):
    
    def __init__(self, matrix_generation_method):
        """
        Class constructor.
        
        @param matrix_generation_method: One of the available matrix generation types available. Currently:
            - load: Load an already created matrix from disk.
            - rmsd: Root Mean Square deviation of one body.
            - distance: Euclidean distance of the geometrical center of one body.
        """
        self.matrix_generation_method = matrix_generation_method
        self.distance_matrix = None
        if not self.matrix_type in ["load","rmsd","distance"]:
            print "[Error] Incorrect matrix creation option: %s"%matrix_generation_method
            exit()
    
    def create_matrix(self, matrix_base_path, trajectory_handler, parameters):
        """
        Generates a matrix with the method used in the handler creation.
        
        @param matrix_base_path: 
        @param parameters: 
        
        @return: The created matrix.
        """
        if self.matrix_generation_method == "load":
            self.distance_matrix = pyRMSD_MatrixHandler.load_matrix(parameters["matrix_path"])
            
        elif self.matrix_generation_method == "rmsd":
            self.distance_matrix =  RMSDMatrixBuilder.build(trajectory_handler)
            
        elif self.matrix_generation_method == "distance":
            self.distance_matrix =  EuclideanDistanceMatrixBuilder().build()
            
        return self.distance_matrix
    
    def save_matrix(self, matrix_save_file):
        """
        Writes matrix contents to disk.
        
        @param matrix_save_file: Complete path (with filename) where to save the matrix.
        """
        self.check_matrix_calculated_error()
        pyRMSD_MatrixHandler.save_matrix(matrix_save_file, self.distance_matrix)
    
    def save_statistics(self, matrix_base_path):
        """
        Writes matrix statistics to disk in JSON format.
        
        @param matrix_base_path: The folder where to save the 'statistics.json' file.
        """
        self.check_matrix_calculated_error()
        pyRMSD_MatrixHandler.save_statistics(matrix_base_path, self.distance_matrix)
    
    def check_matrix_calculated_error(self):
        """
        Exits the program if the matrix wasn't calculated yet.
        """
        if self.distance_matrix is None:
            print "[ERROR][MatrixHandler::save_statistics] Matrix is not been calculated yet."
            exit()
    