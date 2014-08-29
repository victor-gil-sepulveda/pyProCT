"""
Created on 6/8/2014

@author: victor
"""

class DataHandler(object):
    
    def __init__(self, parameters):
        self.data_handler = self.__create_data_handler(parameters)
        self.matrix_handler = self.__calculate_matrix(parameters["matrix"])
    
#     #TODO: cambiar la estructura de los paquetes:
#     data
#         matrix
#             protein
#             array
#         handlers
#             protein
#             array
#     TYPE = family[::special] like protein::single, or protein::ensemble 
#     Handling plugin access: data.handlers.TYPE.SUBTYPE
#     Matrix plugin access: data.matrix.TYPE.plugin+MatrixBuilder
#     Only some handler-matrix are ok.
#     
    def __create_data_handler(self, parameters):
        pass
    
    def __calculate_matrix(self, parameters):
        pass
    
        
    