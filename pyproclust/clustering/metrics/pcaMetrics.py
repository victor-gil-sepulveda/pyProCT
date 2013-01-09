'''
Created on 13/08/2012

@author: victor
'''
import numpy
import pyRMSD.calculators
from pyRMSD import availableCalculators
from pyRMSD.utils.proteinReading import flattenCoords

class PCAMetric(object):
    """
    This is a metric based in the PCA of the clusters composing a clustering,
    """
    def __init__(self):
        pass
    
    def evaluate(self, clustering, trajectory_handler):
        """
        Calculates the value of the PCA metric, which is the mean of the largest eigenvalue obtained from the PCA (the one corresponding
        to the axis of bigger variability) weighted by the cluster size.
        
        @param clustering: The clustering we want to calculate the metric.
        @param trajectory_handler: The TrajectoryHandler containing the coordinates as prody-like coordsets.
        
        @return: the value of the metric.
        """
        # Now do calculation for each one of the clusters
        pca_mean_val = 0.;
        for c in clustering.clusters:
            # Pick the coordinates (ensuring that we are copying them)
            coordinates_of_this_cluster_conformations = numpy.array(trajectory_handler.coordinates[c.all_elements])
            
            # Make an iterative superposition (to get the minimum RMSD of all with respect to a mean conformation)
            coordinates_copy = flattenCoords(coordinates_of_this_cluster_conformations)
            PCAMetric.do_iterative_superposition( coordinates_copy, 
                                                 len(coordinates_of_this_cluster_conformations), 
                                                 trajectory_handler.atoms_per_conformation)
            
            # Calculate the covariance matrix
            coordinates_copy = coordinates_copy.reshape((len(coordinates_of_this_cluster_conformations),trajectory_handler.atoms_per_conformation,3))
            covariance_matrix = PCAMetric.create_covariance_matrix(coordinates_copy)
            
            # And then the eigenvalue we are interested in
            pca_mean_val += PCAMetric.calculate_biggest_eigenvalue(covariance_matrix)
            
        return  pca_mean_val /clustering.total_number_of_elements
    
    @classmethod
    def do_iterative_superposition(cls, flattened_coordinates, number_of_conformations, atoms_per_conformation):
        """
        Performs an iterative superposition over the subset of conformations of a cluster. 
        
        @param flattened_coordinates: The flattened array of these conformations' coordinates. This function
        MODIFIES the array contents.
        @param number_of_conformations: The number of conformations (so the number of elements of the cluster)
        @param atoms_per_conformation: Exactly what the name suggests.
        
        """
        calculator_id = availableCalculators.availableCalculators()["QTRFIT_OMP_CALCULATOR"]
        pyRMSD.calculators.iterativeSuperposition(calculator_id, 
                                                  flattened_coordinates,
                                                  atoms_per_conformation, 
                                                  number_of_conformations,
                                                  8, # OpenMP threads
                                                  0, # CUDA STUFF
                                                  0)
    
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
        for coords in coordinates.reshape((number_of_conformations, coordinates_per_conformation)):
            deviations = coords - mean
            covariance_matrix += numpy.outer(deviations, deviations)
        return covariance_matrix / number_of_conformations
    
    @classmethod
    def calculate_biggest_eigenvalue(cls, covariance_matrix):
        """
        Calculates the eigenvectors and eigenvalues of a covariance matrix.
        
        @param covariance_matrix: 
        """
        values, vectors = numpy.linalg.eigh(covariance_matrix)
        revert = list(range(len(values)-1, -1, -1))
        values = values[revert]
        vectors = vectors[:, revert]
        which = values > 1e-8
        eigvals = values[which]
        print eigvals[0:10]
        return  eigvals[0]