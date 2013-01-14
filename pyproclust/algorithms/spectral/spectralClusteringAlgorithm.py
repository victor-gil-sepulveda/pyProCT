'''
Created on 14/08/2012

@author: victor
'''
from pyproclust.clustering.clustering import Clustering
from pyproclust.clustering.metrics.graphMetrics import d as degree
import numpy
import scipy.linalg
from pyproclust.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from scipy.spatial.distance import pdist
from pyRMSD.condensedMatrix import CondensedMatrix

class SpectralClusteringAlgorithm(object):
    '''
    '''
    
    @staticmethod
    def laplacian_calculation_types(): return  ["PYTHON","NUMPY", "NUMPY_PURE"]
    
    def __init__(self, condensed_matrix, max_clusters, laplacian_calculation_type = "PYTHON"):
        
        print "Calculating W ..."
        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(condensed_matrix)
        
        if not laplacian_calculation_type in SpectralClusteringAlgorithm.laplacian_calculation_types():
            print "[ERROR::SpectralClusteringAlgorithm] Type " ,laplacian_calculation_type, "is not a correct type. Use one of these instead: ", SpectralClusteringAlgorithm.laplacian_calculation_types()
            exit()
        
        L, D = SpectralClusteringAlgorithm.calculate_laplacian(W, condensed_matrix, laplacian_calculation_type)
        
        # Eigenvector i is v[:,i]
        w, vr = scipy.linalg.eig(L, D, right = True, overwrite_a = True, overwrite_b = True)
        
        # Pick the first N rows of the eigenvectors matrix (a row is v[i])
        self.eigen_distances = CondensedMatrix(pdist(vr))
    
    @classmethod
    def calculate_laplacian(cls, W, condensed_matrix, laplacian_calculation_type): 
        if laplacian_calculation_type == "PYTHON":
            print "Calculating D ..."
            D = SpectralClusteringAlgorithm.calculate_degree_matrix(condensed_matrix)
            print "Calculating L ..."
            return SpectralClusteringAlgorithm.calculate_laplacian_python(W, D), numpy.diag(D)
            
        elif laplacian_calculation_type == "NUMPY_PURE":
            print "Calculating D ..."
            D = SpectralClusteringAlgorithm.calculate_degree_matrix(condensed_matrix)
            print "Calculating L ..."
            return SpectralClusteringAlgorithm.calculate_laplacian_numpy_pure(W, D), numpy.diag(D)
            
        elif laplacian_calculation_type == "NUMPY":
            print "Calculating Dinv ..."
            Dinv = numpy.diag(SpectralClusteringAlgorithm.calculate_inverse_degree_matrix(condensed_matrix))
            print "Calculating L ..."
            return SpectralClusteringAlgorithm.calculate_laplacian_numpy(W, Dinv), Dinv
        
        else:
            print "[ERROR SpectralClusteringAlgorithm::calculate_laplacian] Not defined laplacian calculator type: ", laplacian_calculation_type
            exit()
    
    @classmethod
    def calculate_adjacency_matrix(cls, matrix):
        # The adjacency matrix
        data_std = matrix.calculateVariance() 
        
        W_tmp = numpy.zeros((matrix.row_length,)*2, dtype = numpy.float16)
        
        for i in range(matrix.row_length):
            for j in range(i, matrix.row_length):
                W_tmp[i,j] = matrix[i,j]
                W_tmp[j,i] = matrix[i,j]
        
        W = numpy.exp( -W_tmp / data_std).astype(numpy.float16)
        
        # W is squared and symmetric, and its diagonal is 1.
        return W
    
    @classmethod
    def calculate_degree_matrix(cls, matrix):
        # The degree matrix is the diagonal matrix with the degrees for each element.
        D = numpy.zeros(matrix.row_length, dtype = numpy.float16)
        for i in range(matrix.row_length):
            D[i] = degree(i, matrix)
        return D
    
    @classmethod
    def calculate_inverse_degree_matrix(cls, matrix):
        # The degree matrix is the diagonal matrix with the degrees for each element.
        D = numpy.zeros(matrix.row_length, dtype = numpy.float16)
        for i in range(matrix.row_length):
            D[i] = 1./degree(i, matrix)
        return D
    
    @classmethod
    def calculate_laplacian_python(cls, W, D):
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
        # L = I - (D^-1 * W), or  - (D^-1 * W) + I
        L = scipy.linalg.fblas.dgemm(alpha=-1.0, a=Dinv.T, b=W.T, trans_b=True)
        for i in range(L.shape[0]):
            L[i][i] += 1
        return L
    
    @classmethod
    def calculate_laplacian_numpy_pure(cls, W, D):
        I = numpy.identity(W.shape[0]) # Generating I
        Dinv = numpy.linalg.inv(numpy.diag(D).astype(numpy.float32))
        L = I - numpy.dot(Dinv,W)
        return L
            
    def perform_clustering(self, kwargs):
        # Mandatory parameter
        k = int(kwargs["k"])
        
        if k > len(self.eigen_distances): # If k > max_number_of_clusters @ init
            print "[ERROR SpectralClusteringAlgorithm::perform_clustering] this algorithm was defined to generate at maximum ",len(self.eigen_distances)," clusters."
            
        k_medoids_args = {"k":k,
                          "seeding_max_cutoff":-1,
                          "seeding_type": "EQUIDISTANT"}
                          
        k_medoids_alg = KMedoidsAlgorithm(self.eigen_distances)
        
        return k_medoids_alg.perform_clustering(k_medoids_args)
