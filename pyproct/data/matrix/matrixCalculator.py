"""
Created on 2/9/2014

@author: victor
"""
from pyproct.data.matrix.matrixHandler import MatrixHandler
from pyproct.tools.plugins import PluginHandler
import traceback

class MatrixCalculator(object):

    def __init__(self):
        pass
    
    @classmethod
    def calculate(cls, data_handler, matrix_params):
        
        calculator_class = cls.get_calculator(matrix_params)
        
        try:
            distance_matrix = calculator_class.calculate(data_handler, matrix_params)
        except Exception, e:
            print "[ERROR][Driver::postprocess] Impossible to matrix calculation for method: %s"%(calculator.CALCULATION_METHOD)
            print "Message: %s"%str(e)
            traceback.print_exc()
        
        return MatrixHandler(distance_matrix, matrix_params)
    
    @classmethod
    def get_calculator(cls, matrix_params):
        # Get all available calculators
        available_calculators = PluginHandler.get_classes('pyproct.data.matrix', 
                                                          selection_keyword = "MatrixCalculator", 
                                                          skip_list = ["test", "cases"],
                                                          plugin_name = "matrix")
        
        # Choose the calculator we need
        calculation_method = matrix_params["method"]
        calculator_classes = filter(lambda x: x.CALCULATION_METHOD == calculation_method, available_calculators)
        
        if len(calculator_classes) == 0:
            print "[ERROR][MatrixCalculator::calculate] %s is not a registered matrix calculation method."%(calculation_method)
            exit()
        else:
            return calculator_classes[0] 