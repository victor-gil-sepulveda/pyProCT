"""
Created on 1/9/2014

@author: victor
"""
from pyproct.data.matrix.protein.cases.rmsd.cartesiansCase import RMSDMatrixBuilder
from pyproct.data.matrix.protein.cases.rmsd.autoChainMappingCase import ChainMappingBuilder
from pyproct.data.matrix.protein.cases.rmsd.dihedralsCase import DihedralRMSDBuilder

class RMSDMatrixCalculator(object):
    """
    Handles the rmsd for trajectories cases 
    """
    
    CALCULATION_METHOD = "rmsd::ensemble"
    
    def __init__(self, params):
        pass

    @classmethod
    def calculate(cls, data_handler, matrix_params):
        """
        :param parameters: One dictionary entry with at least the keys "method" and
        "parameters":
        
        {
            "method": STRING,
            "parameters":{
                ...
            }
        }
        
        "method": One of the available matrix generation types available. 
        
        For proteins:

        - "rmsd": Root Mean Square deviation of one body
    
                "parameters":{
                    "type": ENUM ["COORDINATES", "DIHEDRALS"]
                    "fit_selection": STRING,
                    "calc_selection": STRING,
                    "calculator_type": ENUM,
                    "calculator_options": OBJECT,
                    "chain_map": BOOL
                }
                
        "type": Type of coordinates used to get the rmsd. It is "COORDINATES" by default.
        "fit_selection": The Prody selection string used to describe the atoms to be superposed.
        "calc_selection": Another Prody selection string that describes the atoms used to calculate RMSD.
        "calculator_type": One of the calculators in pyRMSD.
        "chain_map": Calculates the RMSD of the best mapping of chains e.g. if one has the following 
        tetramer
        
                B
            A       C
                D
    
        and this "reordered" tetramer
    
                A
            B       D
                C
    
        given that all chains are equal, the 'normal' RMSD would be calculated (by default with the 
        chain ordering ABCD Vs ABCD, with a high RMSD value instead of 0, that would be the optimum 
        value for the ordering ABCD Vs BADC.
    
        @return: The created matrix.
        """
        
        coords_type = matrix_params.get_value("type", default_value="COORDINATES")

        if coords_type == "COORDINATES":
            mapping = matrix_params.get_value("chain_map", default_value=False)

            if not mapping:
                return  RMSDMatrixBuilder.build(data_handler, matrix_params)
            else:
                print "Performing Chain Mapping. This may take some time ..."
                return ChainMappingBuilder.calcRMSDMatrix(data_handler.get_data(),
                                matrix_params.get_value("calculator_type", default_value="QCP_SERIAL_CALCULATOR"),
                                matrix_params.get_value("fit_selection", default_value="name CA"))

        elif coords_type == "DIHEDRALS":
            return DihedralRMSDBuilder.build(data_handler.get_data())

