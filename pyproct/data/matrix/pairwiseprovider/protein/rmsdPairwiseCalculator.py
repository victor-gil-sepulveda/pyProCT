'''
Created on Mar 4, 2016

@author: victor
'''

class PairwiseRMSDProvider:
    pass


class PairwiseRMSDCalculator:
    """
    Handles the calculation of pairwise rmsd 
    """
    
    CALCULATION_METHOD = "pairwise::rmsd::ensemble"
    
    def __init__(self, params):
        pass

    @classmethod
    def calculate(cls, data_handler, matrix_params):
        """
        @return: a pairwise rmsd distances provider
        """
        coords_type = matrix_params.get_value("type", default_value="COORDINATES")

        if coords_type == "COORDINATES":
            mapping = matrix_params.get_value("chain_map", default_value=False)

            if not mapping:
                return  PairwiseRMSDProvider(data_handler)
            else:
                print "[ERROR PairwiseRMSDCalculator - Chain Mapping] Not implemented"
                exit()

        elif coords_type == "DIHEDRALS":
            print "[ERROR PairwiseRMSDCalculator - Dihedrals] Not implemented"
            exit()