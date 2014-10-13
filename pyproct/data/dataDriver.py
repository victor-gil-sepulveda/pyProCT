"""
Created on 6/8/2014

@author: victor
"""
import os
from pyproct.data.matrix.matrixCalculator import MatrixCalculator
from pyproct.data.handler.dataHandler import DataHandler
from pyproct.driver.time.timerHandler import timed_method
from pyproct.tools import plotTools

class DataDriver(object):
    
    timer = None
    
    def __init__(self):
        pass
    
    @classmethod
    def run(cls, parameters, workspace_handler, timer, generated_files):
        DataDriver.timer = timer
        
        data_handler = cls.load_data(parameters)
        
        matrix_handler = cls.calc_matrix(data_handler, 
                                         parameters["matrix"])
        
        # Save statistics
        statistics_file_path = matrix_handler.save_statistics(workspace_handler["matrix"])
        generated_files.append({
                                "description":"Matrix statistics",
                                "path":os.path.abspath(statistics_file_path),
                                "type":"text"
        })
        
        # Save matrix contents
        if "filename" in parameters["matrix"]:
            cls.save_matrix(matrix_handler, 
                            workspace_handler,
                            parameters["matrix"])

        # Plot matrix
        if "image" in parameters["matrix"]:
            cls.plot_matrix(matrix_handler, 
                            workspace_handler,
                            parameters["matrix"], 
                            generated_files)
        
        return data_handler, matrix_handler 
    
    @classmethod
    @timed_method(timer, "Data Loading")
    def load_data(cls, parameters):
        return DataHandler(parameters)
        
    @classmethod
    @timed_method(timer, "Matrix Calculation")
    def calc_matrix(cls, data_handler, matrix_parameters):
        return MatrixCalculator.calculate(data_handler, matrix_parameters)
    
    @classmethod
    @timed_method(timer, "Matrix Save")
    def save_matrix(cls, matrix_handler, workspace_handler,  parameters):
        matrix_handler.save_matrix(os.path.join(workspace_handler["matrix"], 
                                                    parameters["filename"]))
    
    @classmethod
    @timed_method(timer, "Matrix Imaging")
    def plot_matrix(cls, matrix_handler, workspace_handler, parameters, generated_files):
        matrix_image_file_path = os.path.join(workspace_handler["matrix"], 
                                              parameters["image"]["filename"])
        max_dim = parameters.get_value("image.dimension", default_value = 1000)
        
        plotTools.matrixToImage(matrix_handler.distance_matrix, 
                                matrix_image_file_path, 
                                max_dim=max_dim)
        
        generated_files.append({
                               "description":"Matrix image",
                               "path":os.path.abspath(matrix_image_file_path),
                               "type":"image"
        })
