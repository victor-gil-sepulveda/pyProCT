"""
Created on 13/02/2013

@author: victor
"""
from pyRMSD.matrixHandler import MatrixHandler as pyRMSD_MatrixHandler
from pyproct.data.proteins.matrix.rmsdMatrixBuilder import RMSDMatrixBuilder
from pyproct.data.proteins.matrix.autoChainMappingRMSDMatrixBuilder import ChainMappingRMSDMatrixCalculator
from pyproct.data.proteins.matrix.dihedralRMSDMatrixCalculator import DihedralRMSDMatrixCalculator
from pyproct.data.proteins.matrix.euclideanDistanceMatrixBuilder import EuclideanDistanceMatrixBuilder

class MatrixHandler(object):

    def __init__(self, matrix_parameters):
        """
        Class constructor.

        @param matrix_parameters

        One JSON entry with this entries:

        {
            matrix: {
                method: method_type,
                parameters:{

                }
            }
        }

        method_type: One of the available matrix generation types available. Currently:

            - 'load': Load an already created matrix from disk

                    parameters:{
                        path: String,
                    }

                    'path' the path where the matrix is going to be loaded.

            - 'rmsd': Root Mean Square deviation of one body

                    parameters:{
                        fit_selection: String,
                        calc_selection: String,
                        calculator_type: Enum
                    }

                    'fit_selection' is the Prody selection string used to describe the atoms to be superposed.
                    'calc_selection' another Prody selection string that describes the atoms used to calculate RMSD.
                    'calculator_type' one of the calculators in pyRMSD.

            - 'auto_chain_map_rmsd': RMSD of the best mapping of chains (e.g. if one has the tetramer
                    B
                A       C
                    D

                and

                    A
                B       D
                    C

                and all chains are equal, the 'normal' RMSD would be calculated with ABCD vs BADC giving place to a high RMSD value
                instead of 0, that would be the optimum value.

                    parameters:{
                        fit_selection: String,
                        calculator_type: Enum
                    }

                    'fit_selection' is the Prody selection string used to describe the atoms to be superposed.
                    'calculator_type' one of the calculators in pyRMSD.

            - 'distance': Euclidean distance of the geometrical center of one body.

                    parameters:{
                        fit_selection:  String,
                        body_selection: String,
                    }

                    'fit_selection' is the Prody selection string used to describe the atoms to be superposed.
                    'body_selection' another Prody selection string that describes the element that will be used
                    to get the euclidean distances.
        """
        self.matrix_parameters = matrix_parameters

        self.distance_matrix = None

        if not self.matrix_parameters["method"] in ["load","rmsd","distance"]:
            print "[Error] Incorrect matrix creation option: %s"%self.matrix_parameters["method"]
            exit()

    def calculate_matrix(self, trajectory_handler):
        """
        Calculates a matrix with the method used in the handler creation.

        @param matrix_base_path:
        @param parameters:

        @return: The created matrix.
        """
        if self.matrix_parameters["method"] == "load":
            self.distance_matrix = pyRMSD_MatrixHandler.load_matrix(self.matrix_parameters["parameters"]["path"])

        elif self.matrix_parameters["method"] == "rmsd":
            type = self.matrix_parameters.get_value("parameters.type", default_value="COORDINATES")

            if type == "COORDINATES":
                mapping = self.matrix_parameters.get_value("parameters.chain_map", default_value=False)

                if not mapping:
                    self.distance_matrix =  RMSDMatrixBuilder.build(trajectory_handler, self.matrix_parameters["parameters"])
                else:
                    print "Performing Chain Mapping. This may take some time ..."
                    self.distance_matrix =  ChainMappingRMSDMatrixCalculator.calcRMSDMatrix(trajectory_handler.getMergedStructure(),
                                    self.matrix_parameters.get_value("parameters.calculator_type", default_value="QCP_SERIAL_CALCULATOR"),
                                    self.matrix_parameters.get_value("parameters.fit_selection", default_value="name CA"))

            elif type == "DIHEDRALS":
                self.distance_matrix = DihedralRMSDMatrixCalculator.build(trajectory_handler.getMergedStructure())
                pass
            else:
                print self.matrix_parameters
                print "[Error] Incorrect matrix 'rmsd' creation option: %s"%type
                exit()
        elif self.matrix_parameters["method"] == "distance":
            self.distance_matrix =  EuclideanDistanceMatrixBuilder.build(trajectory_handler,self.matrix_parameters["parameters"])

        return self.distance_matrix

    def save_matrix(self, matrix_path):
        """
        Writes matrix contents to disk.

        @param matrix_save_file: Complete path (with filename) where to save the matrix.
        """
        self.check_matrix_calculated_error()
        pyRMSD_MatrixHandler.save_matrix(matrix_path, self.distance_matrix)

    def save_statistics(self, matrix_base_path):
        """
        Writes matrix statistics to disk in JSON format.

        @param matrix_base_path: The folder where to save the 'statistics.json' file.
        """
        self.check_matrix_calculated_error()
        return pyRMSD_MatrixHandler.save_statistics(matrix_base_path, self.distance_matrix)

    def check_matrix_calculated_error(self):
        """
        Exits the program if the matrix wasn't calculated yet.
        """
        if self.distance_matrix is None:
            print "[ERROR][MatrixHandler::save_statistics] Matrix has not been calculated yet."
            exit()

