'''
Created on 13/08/2012

@author: victor
'''
import numpy
from pyRMSD.RMSDCalculator import RMSDCalculator
from pyproclust.driver.handlers.timerHandler import TimerHandler
import scipy.linalg

class PCAMetric(object):
    """
    This is a metric based in the PCA of the clusters composing a clustering,
    """
    def __init__(self, trajectory_handler):
        """
        @param trajectory_handler: The TrajectoryHandler containing the coordinates as prody-like coordsets.
        """
        self.copied_coordinates = trajectory_handler.getWorkingCoordinates()

    def evaluate(self, clustering):
        """
        Calculates the value of the PCA metric, which is the mean of the largest eigenvalue obtained from the PCA (the one corresponding
        to the axis of bigger variability) weighted by the cluster size.
        
        @param clustering: The clustering we want to calculate the metric.
        
        @return: the value of the metric.
        """
        # Pca for each one of the clusters
        pca_mean_val = 0.;
        
        for c in clustering.clusters:
            # Pick the coordinates (ensuring that we are copying them)
            coordinates_of_this_cluster_conformations = self.copied_coordinates[c.all_elements]
            calculator = RMSDCalculator(coordsets = coordinates_of_this_cluster_conformations,
                                        calculatorType = "QTRFIT_OMP_CALCULATOR")
            
            # Make an iterative superposition (to get the minimum RMSD of all with respect to a mean conformation)
            calculator.iterativeSuperposition()

            # Calculate the covariance matrix
            covariance_matrix = PCAMetric.create_covariance_matrix(coordinates_of_this_cluster_conformations)
            
            # And then the eigenvalue we are interested in
            pca_mean_val += PCAMetric.calculate_biggest_eigenvalue(covariance_matrix)

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