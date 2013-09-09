'''
Created on 14/08/2012

@author: victor
'''
import numpy
import scipy.linalg
from pyproclust.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from scipy.spatial.distance import pdist
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproclust.clustering.cluster import gen_clusters_from_class_list
import scipy.cluster.vq
from pyproclust.clustering.clustering import Clustering

class SpectralClusteringAlgorithm(object):
    '''
    Implementation of Normalized Spectral clustering with Lrw Laplacian (Shi and Malik 2000).
    It tries to be both fast and memory efficient in order to be able to work with very big datasets.
    '''
    
    @staticmethod
    def laplacian_calculation_types():
        """
        Returns a list with the correct types of W calculation.
        - PYTHON calculator is fast for large datasets and is the one with lower memory needs.
        - NUMPY Does the calculations using numpy to get some speedup. Is the faster one for small datasets.
        - NUMPY pure uses the real algebraic representation of the calculations. It depends only on Numpy. Useful as golden data
        generator.
        """
        return  ["PYTHON","NUMPY", "NUMPY_PURE"]
    
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
        
        try:
            verbose = kwargs["verbose"]
        except KeyError:
            verbose = False

        if verbose: print "Creating spectral."

        try:
            self.max_clusters = kwargs["max_clusters"]
        except KeyError:
            self.max_clusters = condensed_matrix.row_length
        
        
        try:
            self.sigma_sq = kwargs["sigma_sq"]
            W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(condensed_matrix, self.sigma_sq)
        except KeyError:
            W, self.sigma_sq = self.do_sigma_estimation(condensed_matrix)
            if verbose: print "Sigma estimation: ", self.sigma_sq
        
        try:
            store_W = kwargs["store_W"]
        except KeyError:
            store_W = False
        

        try:
            laplacian_calculation_type = kwargs["laplacian_calculation_type"]
            if not laplacian_calculation_type in SpectralClusteringAlgorithm.laplacian_calculation_types():
                print "[ERROR::SpectralClusteringAlgorithm] Type " ,laplacian_calculation_type, "is not a correct type. Use one of these instead: ", SpectralClusteringAlgorithm.laplacian_calculation_types()
                exit()
        except KeyError:
            laplacian_calculation_type = "NUMPY_PURE"
                
        if verbose: print "Calculating W ..."
        
        # Zero all negative values (similarities cannot be < 0 )
        W[W<0] = 0.0        
        
        if store_W:
            self.W = numpy.copy(W)
        
        L, D = SpectralClusteringAlgorithm.calculate_laplacian(W, condensed_matrix, laplacian_calculation_type, verbose)
        
        # Eigenvector i is v[:,i]
        eigenvalues, self.eigenvectors = scipy.linalg.eig(L, D, right = True, overwrite_a = True, overwrite_b = True)
        # We can try with scipy.sparse.linalg.eigs as the matrix is sparse
#         eigenvalues, self.eigenvectors = scipy.sparse.linalg.eigs(L, k = self.max_clusters, M = D,which ='SM',return_eigenvectors = True)
        
        # Order eigenvectors by eigenvalue (from lowest to biggest)
        idx = eigenvalues.real.argsort()
        self.eigenvectors = self.eigenvectors[:, idx]
        
        # We'll only store the vectors we need, usually << N
        self.eigenvectors = self.eigenvectors[:,:self.max_clusters]
        if verbose: print "Spectral finished."
            
    def perform_clustering(self, kwargs):
        """
        Does the actual clustering by doing a k-medoids clustering of the first k eigenvector rows.
        
        @param kwargs: Dictionary with this mandatory keys:
            - 'k': Number of clusters to generate. Must be <= than max_clusters
            
        @return: a Clustering instance with the clustered data.
        """
        # Mandatory parameter
        k = int(kwargs["k"])
        
        try:
            use_k_medoids = bool(kwargs["use_k_medoids"])
        except KeyError:
            use_k_medoids = True
        
        if k > self.max_clusters: # If k > max_number_of_clusters @ init
            print "[ERROR SpectralClusteringAlgorithm::perform_clustering] this algorithm was defined to generate at maximum ",self.max_clusters," clusters."

        algorithm_details = "Spectral algorithm with k = %d and sigma squared = %.3f" %(int(k), self.sigma_sq)
        
        if use_k_medoids:
            # The row vectors we have are in R^k (so k length)
            eigen_distances = CondensedMatrix(pdist(self.eigenvectors[:,:k]))
            k_medoids_args = {"k":k,
                              "seeding_max_cutoff":-1,
                              "seeding_type": "EQUIDISTANT"}
            k_medoids_alg = KMedoidsAlgorithm(eigen_distances)
            clustering = k_medoids_alg.perform_clustering(k_medoids_args)
            clustering.details = algorithm_details
            return k_medoids_alg.perform_clustering(k_medoids_args)
        else:
            centroid, labels = scipy.cluster.vq.kmeans2(self.eigenvectors[:,:k], 
                                                        k, iter = 1000, minit = 'points')
            del centroid
            clusters = gen_clusters_from_class_list(labels)
            return Clustering(clusters,details = algorithm_details)
    
    @classmethod
    def calculate_laplacian(cls, W, condensed_matrix, laplacian_calculation_type, verbose = False):
        """
        Calculates the RW laplacian: L = I - D^-1 * W 
        
        @param W: Adjacency matrix.
        @param condensed_matrix: The distance matrix for this dataset.
        @param laplacian_calculation_type: One of the calculation types returned by laplacian_calculation_types()
        
        @return: The laplacian.
        """
        if laplacian_calculation_type == "PYTHON":
            if verbose: print "Calculating D ..."
            D = SpectralClusteringAlgorithm.calculate_degree_matrix(W)
            if verbose: print "Calculating L ..."
            return SpectralClusteringAlgorithm.calculate_laplacian_python(W, D), numpy.diag(D)
            
        elif laplacian_calculation_type == "NUMPY":
            if verbose: print "Calculating Dinv ..."
            Dinv = numpy.diag(SpectralClusteringAlgorithm.calculate_inverse_degree_matrix(W))
            if verbose: print "Calculating L ..."
            return SpectralClusteringAlgorithm.calculate_laplacian_numpy(W, Dinv), Dinv
        
        elif laplacian_calculation_type == "NUMPY_PURE":
            if verbose: print "Calculating D ..."
            D = SpectralClusteringAlgorithm.calculate_degree_matrix(W)
            if verbose: print "Calculating L ..."
            return SpectralClusteringAlgorithm.calculate_laplacian_numpy_pure(W, D), numpy.diag(D)
            
        else:
            print "[ERROR SpectralClusteringAlgorithm::calculate_laplacian] Not defined laplacian calculator type: ", laplacian_calculation_type
            exit()
    
    def do_sigma_estimation(self,matrix):
        """
        Does a first sigma approximation so that at most a 10% of the elements are lower than 0.
        
        @param matrix: The distance matrix for this dataset.
        
        @return: The adjacency matrix with the chosen sigma estimation.
        
        """
        sq_first_approx = matrix.calculateVariance() / matrix.calculateMax() # Possible lower bound
        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix( matrix, sq_first_approx)
        number_of_elems = float(matrix.row_length ** 2)
        number_of_negative_elements = len(W[W<0].flatten())
        print "Negative percent",number_of_negative_elements,number_of_elems,number_of_negative_elements / number_of_elems 
        while number_of_negative_elements / number_of_elems > 0.1 :
            sq_first_approx += 0.5
            W = SpectralClusteringAlgorithm.calculate_adjacency_matrix( matrix, sq_first_approx)
            number_of_negative_elements  = float(len(W[W<0].flatten()))
            print "Negative percent",number_of_negative_elements/ number_of_elems 
        
        return W, sq_first_approx
        
    @classmethod
    def calculate_adjacency_matrix(cls, matrix, sigma_sq):
        """
        Calculates the adjacency matrix (W) representing the distance graph for this dataset. It constructs a 
        'fully connected graph', as explained in (Luxburg 2007).
        
        @param matrix: The distance matrix for this dataset.
        
        @param sigma_sq: Is a cutoff parameter for edge creation.
        
        @return: The adjacency matrix.
        """
        # The adjacency matrix
        W_tmp = numpy.zeros((matrix.row_length,)*2, dtype = numpy.float16)
        
        for i in range(matrix.row_length):
            for j in range(i, matrix.row_length):
                W_tmp[i,j] = matrix[i,j]
                W_tmp[j,i] = matrix[i,j]

        W = numpy.exp( - ((W_tmp**2) / float(sigma_sq))).astype(numpy.float16)
        
        # W is squared and symmetric, and its diagonal is 1.
        return W
    
    @classmethod
    def calculate_degree_matrix(cls, W):
        """
        Calculates... the degree matrix (as an array containing the diagonal).
        
        @param W: The adjacency matrix.
        
        @return: The degree matrix as an array.
        """
        # The degree matrix is the diagonal matrix with the degrees for each element.
        # The degree is the sum of the row of W (or column as it's symmetric) for each element.
        return [sum(Wi) for Wi in W]
    
    @classmethod
    def calculate_inverse_degree_matrix(cls, W):
        """
        Calculates the inverse of the degree matrix.
        @see: calculate_degree_matrix
        
        @param W: The adjacency matrix.
        
        @return: The inverse degree matrix as an array.
        """
        return [1./sum(Wi) for Wi in W]
    
    @classmethod
    def calculate_laplacian_python(cls, W, D):
        """
        Implementation of the laplacian calculation for the 'PYTHON' strategy.
        
        @param W: Adjacency matrix.
        @param D: Degree matrix.
        
        @return: Laplacian matrix.
        """
        # L = I - D^-1 * W 
        # Run trough the rows and divide row i by D[i]. Because of construction, diagonal is 1 - (1*(1/D[i]))
        for i in range(W.shape[0]):
            for j in range(W.shape[1]):
                if(i==j):
                    W[i][i] = 1. - (1. / D[i])
                else:
                    W[i][j] = - (W[i][j] / D[i])
        return W
    
    @classmethod
    def calculate_laplacian_numpy(cls, W, Dinv):
        """
        Implementation of the laplacian calculation for the 'NUMPY' strategy.
        
        @param W: Adjacency matrix.
        @param D: Degree matrix.
        
        @return: Laplacian matrix.
        """
        # L = I - (D^-1 * W), or  - (D^-1 * W) + I
        L = scipy.linalg.fblas.dgemm(alpha=-1.0, a=Dinv.T, b=W.T, trans_b=True)
        for i in range(L.shape[0]):
            L[i][i] += 1
        return L
    
    @classmethod
    def calculate_laplacian_numpy_pure(cls, W, D):
        """
        Implementation of the laplacian calculation for the 'NUMPY_PURE' strategy.
        
        @param W: Adjacency matrix.
        @param D: Degree matrix.
        
        @return: Laplacian matrix.
        """
        I = numpy.identity(W.shape[0]) # Generating I
        Dinv = numpy.linalg.inv(numpy.diag(D).astype(numpy.float32))
        L = I - numpy.dot(Dinv,W)
        return L
            
