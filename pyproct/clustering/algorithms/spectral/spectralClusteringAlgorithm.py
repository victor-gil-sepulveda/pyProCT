"""
Created on 14/08/2012

@author: victor
"""
import numpy
from scipy.spatial.distance import pdist
from pyproct.clustering.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
import scipy.cluster.vq
from pyproct.clustering.clustering import Clustering
from pyproct.clustering.cluster import gen_clusters_from_class_list
from pyRMSD.condensedMatrix import CondensedMatrix
import pyproct.clustering.algorithms.spectral.cython.spectralTools as SpectralTools

#
# first -> similarity graph -> fully connected graph (undirected)  -> similarity function exp(- d^2 /2sigma^2)
# second -> weighted adjacency matrix comes from similarity graph (wij = sij)
# first k eigenvectors == k smallest eigenvectors
# Unnorm L= D-W
# L and Lsym are symmetric and positive definite, we can force sparsity (as neigbour versions should have)


# TODO: Eigenvalues 0 are a good indicator of the number of connected components, however to use it, one has to precalculate L. This is not
# efficient with the current architecture. <- PARAMETERS


class SpectralClusteringAlgorithm(object):
    """
    Implementation of Normalized Spectral clustering with Lrw Laplacian (Shi and Malik 2000).
    It tries to be both fast and memory efficient in order to be able to work with very big datasets.
    """
    spectral_types_enum =  ["UNNORMALIZED","NORMALIZED"]



    def __init__(self, condensed_matrix, **kwargs):
        """
        Constructor. Calculates the eigenvectors given a dataset distance matrix. The eigenvector distances would be the
        common part for clusterings with different k.

        @param condensed_matrix: The distance matrix of the dataset.
        @param sigma_sq: The squared value of sigma for the adjacency matrix calculation. If None, the value will be automatically
        calculated.
        @param max_clusters: Maximum number of clusters we will try with this algorithm (for instance with max_clusters = 10
        we can try with ks in range [1..10]
        @param laplacian_calculation_type: The type of calculation.
        @param store_W: If True the object stores the adjacency matrix. Useful for testing.
        @param verbose: If True some messages will be printed.
        """
        self.handle_params(kwargs, max_clusters_default = condensed_matrix.row_length-1)

        print "Initializing Spectral clustering. This may take some time ..."
#         self.verbose = True

        if self.sigma_sq is not None:
            if self.verbose: print "Calculating W with sigma = %f estimation..."%self.sigma_sq
            W = SpectralTools.calculate_fully_connected_adjacency_matrix(condensed_matrix, self.sigma_sq)
        else:
            if self.verbose: print "Calculating W with sigma estimation..."
            sigmas = SpectralTools.local_sigma_estimation(condensed_matrix)
            W  = SpectralTools.calculate_fully_connected_adjacency_matrix_with_sigma_estimation(condensed_matrix, sigmas)
            self.sigma_sq = numpy.mean(sigmas)**2

        if self.verbose: print "Sigma^2 estimation (mean of local sigmas): ", self.sigma_sq

        if self.force_sparse:
            SpectralTools.force_sparsity(W)

        if self.store_W:
            if self.verbose: print "Storing W ..."
            self.W = numpy.copy(W)

        if self.verbose: print "Calculating Degree Matrix ..."
        D = SpectralTools.calculate_degree_matrix(W)

        if self.verbose: print "Calculating Laplacian ..."
        L =  SpectralTools.calculateUnnormalizedLaplacian(W, D)

        if self.verbose: print "Calculating Eigenvectors ..."
        if self.spectral_type == "UNNORMALIZED":
            v = SpectralTools.calculateUnnormalizedEigenvectors(L, self.max_clusters, self.force_sparse)
        elif self.spectral_type == "NORMALIZED":
            v = SpectralTools.calculateNormalizedEigenvectors(L, D, self.max_clusters, self.force_sparse)

        self.eigenvectors = v # eigenvectors in columns. We need the rows of this matrix for the clustering.
        if self.verbose: print "Spectral initialization finished."

    def perform_clustering(self, kwargs):
        """
        Does the actual clustering by doing a k-medoids clustering of the first k eigenvector rows.

        @param kwargs: Dictionary with this mandatory keys:
            - 'k': Number of clusters to generate. Must be <= than max_clusters

        @return: a Clustering instance with the clustered data.
        """
        # Mandatory parameter
        k = int(kwargs["k"])

        if k > self.max_clusters:
            print "[ERROR SpectralClusteringAlgorithm::perform_clustering] this algorithm was defined to generate at most %d clusters."%self.max_clusters,

        algorithm_details = "Spectral algorithm with k = %d and sigma squared = %.3f" %(int(k), self.sigma_sq)

        if self.use_k_medoids:
            # The row vectors we have are in R^k (so k length)
            eigen_distances = CondensedMatrix(pdist(self.eigenvectors[:,:k]))
            k_medoids_args = {
                              "k":k,
                              "seeding_max_cutoff":-1,
                              "seeding_type": "RANDOM"
                              }

            k_medoids_alg = KMedoidsAlgorithm(eigen_distances)
            clustering = k_medoids_alg.perform_clustering(k_medoids_args)
            clustering.details = algorithm_details
            return k_medoids_alg.perform_clustering(k_medoids_args)
        else:
            centroid, labels = scipy.cluster.vq.kmeans2(self.eigenvectors[:,:k],
                                                        k,
                                                        iter = 1000,
                                                        minit = 'random')
            del centroid
            clusters = gen_clusters_from_class_list(labels)
            return Clustering(clusters,details = algorithm_details)

    def handle_params(self, params, max_clusters_default):
        try:
            self.verbose = params["verbose"]
        except KeyError:
            self.verbose = False

        try:
            self.max_clusters = int(params["max_clusters"])
        except KeyError:
            self.max_clusters = max_clusters_default

        try:
            self.sigma_sq = float(params["sigma_sq"])
            self.sigma_estimation = False
        except KeyError:
            self.sigma_sq = None
            self.sigma_estimation = True

        try:
            self.store_W = params["store_W"]
        except KeyError:
            self.store_W = False

        try:
            self.spectral_type = params["type"]
            if not self.spectral_type in SpectralClusteringAlgorithm.spectral_types_enum:
                print "[ERROR::SpectralClusteringAlgorithm] Type " ,self.spectral_type,\
                "is not a correct type. Use one of these instead: ", SpectralClusteringAlgorithm.spectral_types_enum
                exit()
        except KeyError:
            self.spectral_type = "UNNORMALIZED"

        try:
            self.force_sparse = params["force_sparse"]
        except KeyError:
            self.force_sparse = False

        try:
            self.use_k_medoids = bool(params["use_k_medoids"])
        except KeyError:
            self.use_k_medoids = True

