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
        # Get all available calculators
        available_calculators = PluginHandler.get_classes('pyproct.data.matrix', 
                                                          selection_keyword = "MatrixCalculator", 
                                                          skip_list = ["test", "cases"],
                                                          plugin_name = "matrix")
        
        # Choose the calculator we need
        calculator_class = None
        calculation_method = matrix_params["method"]
        for calculator in available_calculators:
            if calculator.CALCULATION_METHOD == calculation_method:
                calculator_class = calculator  
        
        if calculator_class is None:
            print "[ERROR][MatrixCalculator::calculate] %s is not a registered matrix calculation method."%(calculation_method)
            exit()
            
        try:
            distance_matrix = calculator_class.calculate(data_handler, matrix_params)
        except Exception, e:
            print "[ERROR][Driver::postprocess] Impossible to perform '%s' postprocessing action."%(calculator.CALCULATION_METHOD)
            print "Message: %s"%str(e)
            traceback.print_exc()
        
        return MatrixHandler(distance_matrix, matrix_params)