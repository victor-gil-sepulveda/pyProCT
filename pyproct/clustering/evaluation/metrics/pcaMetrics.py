"""
Created on 13/08/2012

@author: victor
"""
import numpy
from pyRMSD.RMSDCalculator import RMSDCalculator
from pyproct.driver.time.timerHandler import TimerHandler
import scipy.linalg

class PCAMetric(object):
    """
    This is a metric based in the PCA of the clusters composing a clustering,
    """
    def __init__(self, trajectory_handler):
        """
        @param trajectory_handler: The TrajectoryHandler containing the coordinates as prody-like coordsets.
        """
        self.fitting_coordinates = trajectory_handler.getFittingCoordinates()
        self.calculation_coordinates = trajectory_handler.getCalculationCoordinates()

    def evaluate(self, clustering):
        """
        Calculates the value of the PCA metric, which is the mean of the largest eigenvalue obtained from the PCA (the one corresponding
        to the axis of bigger variability) weighted by the cluster size.
        
        @param clustering: The clustering we want to calculate the metric.
        
        @return: the value of the metric.
        """
        # Pca for each one of the clusters
        pca_mean_val = 0.;
        MAX_ELEMENTS = 1000
        for c in clustering.clusters:
            # Pick the coordinates (ensuring that we are copying them)
            element_indexes = c.all_elements
            ###################
            # Performance hack
            ###################
            # As it can be very slow for big clusters (i.e. > 3k elements) we'll compress this clusters 
            # before calculating PCA. It should increase variance but will allow calculations.
            # It should use the kmedoids compressor
            if len(c.all_elements) > MAX_ELEMENTS:
                element_indexes = c.get_random_sample(MAX_ELEMENTS)
                print "[PCA] Random sampling too big cluster to improve performance (%d elements -> %d elements)."%(len(c.all_elements),MAX_ELEMENTS)
            ###################
            
            fitting_coordinates_of_this_cluster = self.fitting_coordinates[element_indexes]
            
            calculator = RMSDCalculator(calculatorType = "QTRFIT_SERIAL_CALCULATOR",
                                        fittingCoordsets = fitting_coordinates_of_this_cluster)
            
            if self.calculation_coordinates is not None:
                calculation_coordinates_of_this_cluster = self.calculation_coordinates[element_indexes]
                calculator = RMSDCalculator(calculatorType = "QTRFIT_SERIAL_CALCULATOR",
                                            fittingCoordsets = fitting_coordinates_of_this_cluster,
                                            calculationCoordsets = calculation_coordinates_of_this_cluster)
            
            # Make an iterative superposition (to get the minimum RMSD of all with respect to a mean conformation)
            calculator.iterativeSuperposition()

            # Calculate the covariance matrix
            if self.calculation_coordinates is None:
                covariance_matrix = PCAMetric.create_covariance_matrix(fitting_coordinates_of_this_cluster)
            else:
                covariance_matrix = PCAMetric.create_covariance_matrix(calculation_coordinates_of_this_cluster)
            
            # And then the eigenvalue we are interested in
            pca_mean_val += PCAMetric.calculate_biggest_eigenvalue(covariance_matrix)
        print "PCA finished"
        return  pca_mean_val /clustering.total_number_of_elements
    
    @classmethod
    def create_covariance_matrix(cls, coordinates):
        """
        Calculates the covariance matrix for a given number of conformations
        
        @param coordinates: Coordinates of the trajectory frames we want to use in the covariance matrix calculation.
        This coordinates are stored in the format [conformation 1..N] where conformation is [Atom 1...N] Where atom is [x,y,z],
        so it's a (C,A,3) shaped matrix where C is the number of conformations and N the number of atoms per conformation.

        @return: The covariance matrix.
        """
        number_of_conformations = coordinates.shape[0]
        number_of_atoms =  coordinates.shape[1]
        coordinates_per_conformation = number_of_atoms * 3
        covariance_matrix = numpy.zeros((coordinates_per_conformation, coordinates_per_conformation))
        coordinates = coordinates.reshape((number_of_conformations, coordinates_per_conformation))
        # Mean structure
        mean = coordinates.mean(0)
        # Changed for efficiency
        for coords in coordinates:
            deviations = coords - mean
            covariance_matrix += numpy.outer(deviations, deviations)
        return covariance_matrix / number_of_conformations
    
    @classmethod
    def calculate_biggest_eigenvalue(cls, covariance_matrix):
        """
        Calculates the eigenvectors and eigenvalues of a covariance matrix.
        
        @param covariance_matrix: The covariance matrix.
        
        @return: The first (bigger) eigenvalue, which gives an idea of the relative importance of the 
        axis responsible of most of the variance. 
        """
        timer = TimerHandler()
        timer.start("eigen2")
        eigvals = scipy.linalg.eigh(covariance_matrix, 
                                   eigvals_only = True, 
                                   eigvals = (covariance_matrix.shape[0] -1,covariance_matrix.shape[0]-1), 
                                   overwrite_a = True)
        return  eigvals[0]
#         THIS IS THE MOST COSTLY OPERATION
#         values, vectors = numpy.linalg.eigh(covariance_matrix)
#         timer.stop("eigen2")
#         print "INSIDE EIGENVALUES"
#         print timer
#         revert = list(range(len(values)-1, -1, -1))
#         values = values[revert]
#         vectors = vectors[:, revert]
#         which = values > 1e-8
#         eigvals = values[which]
#         print eigvals
#         return  eigvals[0]