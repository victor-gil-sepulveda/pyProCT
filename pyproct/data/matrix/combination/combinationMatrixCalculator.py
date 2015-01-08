"""
Created on 27/11/2014
 
@author: victor
"""
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.data.matrix.matrixCalculator import MatrixCalculator
 
class combinationMatrixCalculator(object):
    """
    Handles the rmsd for trajectories cases 
    """
     
    CALCULATION_METHOD = "matrix::combination"
     
    def __init__(self, params):
        pass
 
    @classmethod
    def calculate(cls, data_handler, matrix_params):
        """
        :param parameters: a ProtocolParameters object with at least one "matrices" 
        and "combination" attributes. 
         
        Inside the "matrices" object any type of "matrix" object can be defined. The 
        attribute id would be its id in the "combination" clause.
         
        Ex:
        "matrix":{
            "method":"matrix::combination",
            "parameters":{
                "matrices":{
                    "loaded_matrix":{
                        "method": "matrix::load",
                        "parameters": {
                            ...
                        }
                    },
                    "rmsd_matrix":{
                        "method": "rmsd:ensemble",
                        "parameters": {
                            "calculator_type": "QCP_OMP_CALCULATOR",
                            "fit_selection": "chain A",
                            ...
                        }
                    }
                },
                "combination": ["add",[]]
            }
         
        }
     
        @return: The created matrix.
        """
        forbidden = ["add","sub","mult"]
        matrices_descr = matrix_params.get_value("matrices", default_value={})
        matrices = {}
        operations = matrix_params.get_value("combination", default_value="")
         
        number_of_matrix_ids = len(matrices_descr.keys()) 
        if number_of_matrix_ids == 0:
            #raise 
            pass
        else:
            #check that there are no forbidden ids
            for matrix_id in matrices_descr:
                if matrix_id in forbidden:
                    raise KeyError("Forbidden matrix id.")
         
        # Load the matrices
        for matrix_id in matrices_descr:
            print "Calculating %s"%(matrix_id)
            matrices[matrix_id] = MatrixCalculator.calculate(data_handler, 
                                                             matrices_descr[matrix_id])
             
        # Do the combination
        return combine(operations, matrices)
 
 
def combine(operation, matrices):
    try:
        operator = operation[0]
    except:
        operator = operation
     
    if operator == "add":
        return add_matrices(combine(operation[1], matrices), 
                            combine(operation[2], matrices))
         
    elif operator == "sub":
        return sub_matrices(combine(operation[1], matrices), 
                            combine(operation[2], matrices))
          
    elif operator == "mult":
        return multiply_by_scalar(combine(operation[1], matrices), 
                            combine(operation[2], matrices))
         
    elif operator in matrices:
        # Then is a matrix id
        return matrices[operator]
     
    else:
        # Try to convert it to a float
        try:
            return float(operator)
        except:
            print "[Error combinationMatrixCalculator:combine] Unexpected operator or id (%s)."%operator
            exit()
 
def add_matrices(matrix1, matrix2):
    return CondensedMatrix(matrix1.get_data()+matrix2.get_data())
 
def sub_matrices(matrix1, matrix2):
    return CondensedMatrix(matrix1.get_data()-matrix2.get_data())
 
def multiply_by_scalar(op1, op2):
    is_scalar = [isinstance(op1, float), isinstance(op2, float)]
     
    if is_scalar[0] and not is_scalar[1]:
        scalar = op1
        matrix = op2
    elif not is_scalar[0] and  is_scalar[1]:
        scalar = op2
        matrix = op1
    else:
        print "[Error combinationMatrixCalculator:multiply_by_scalar] One (and only one) operand must be ."
        exit()
         
    return CondensedMatrix(scalar*matrix.get_data())
         
             