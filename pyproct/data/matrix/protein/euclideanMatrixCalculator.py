"""
Created on 13/02/2013

@author: victor
"""
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.data.matrix.protein.cases.euclidean.dihedralsCase import DihedralEuclideanDistanceBuilder
from pyproct.data.matrix.protein.cases.euclidean.cartesiansCase import euclideanDistanceBuilder

class EuclideanMatrixCalculator(object):
    
    CALCULATION_METHOD = "euclidean_distance::ensemble"
    
    def __init__(self):
        pass

    @classmethod
    def calculate(cls, data_handler, matrix_parameters):
        """
        Will generate the CondensedMatrix filled with the all vs all geometric center distances of the "body_selection"
        coordinates (which will usually be a ligand).

        @param trajectory_handler: The handler containing selection strings, pdb info and coordsets.
        @param matrix_parameters: The creation parameters (from the initial script).

        @return: The created distance matrix.
        """
        
        coords_type = matrix_parameters.get_value("type", default_value="COORDINATES")
        builder = None
        
        if coords_type == "COORDINATES":
            print "using coords"
            builder = euclideanDistanceBuilder
            
        elif coords_type == "DIHEDRALS":
            print "using dihedrals"
            builder = DihedralEuclideanDistanceBuilder

        coordinates = builder.build(data_handler,  matrix_parameters)
        distances = builder.calc_distances(coordinates)
        
        return CondensedMatrix(distances)

    
