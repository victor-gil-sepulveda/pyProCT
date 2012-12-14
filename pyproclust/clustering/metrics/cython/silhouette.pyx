# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
'''
Created on 12/07/2012

@author: victor
'''

import numpy
cimport numpy
cimport cython

cdef inline double double_max(double a, double b): return a if a >= b else b
cdef inline double double_min(double a, double b): return a if a <= b else b
DTYPE = numpy.double
ctypedef numpy.double_t DTYPE_t

cdef class CythonSilhouetteCoefficientCalculator(object):
    def __init__(self): 
        pass
    
    @cython.boundscheck(False)
    def evaluate(self, clustering, condensed_distance_matrix):
        cdef int total_number_of_elements,i 
        cdef numpy.ndarray[DTYPE_t, ndim=1] cluster_silhouettes
        
        # Undetermined case
        if (len(clustering.clusters) == 1):
            return 0.
        
        total_number_of_elements = clustering.total_number_of_elements
        
        cluster_silhouettes = numpy.empty(total_number_of_elements)
        
        i = 0
        for cluster in clustering.clusters:
            for element in cluster.all_elements:
                cluster_silhouettes[i] = self.__one_element_silhouette(element,cluster,clustering,condensed_distance_matrix)
                i = i + 1
        
        return numpy.mean(cluster_silhouettes)
    
    @cython.boundscheck(False)
    def __one_element_silhouette(self,int element,cluster,clusterization,condensed_distance_matrix):
        cdef double a_i, b_i
        cdef int where_am_i, i 
        
        a_i = self.__get_average_distance_with_my_cluster(element,cluster,condensed_distance_matrix)
        #where_am_i = clusterization.cluster_index(cluster)
        #################
        # TOTAL REWRITE
        #################
        where_am_i = 0
        for i in range(len(clusterization.clusters)):
            if element in clusterization.clusters[i].all_elements:
                where_am_i = i
                break
        ####### 
        b_i = numpy.finfo(numpy.float32).max
        for i in range(len(clusterization.clusters)):
            
            if where_am_i != i:
                b_i = double_min(b_i,self.__get_average_distance_with_cluster(element,clusterization.clusters[i],condensed_distance_matrix))
        return (b_i-a_i)/double_max(a_i,b_i)
    
    @cython.boundscheck(False)
    cpdef double __sum_cluster_distances(self,int element,cluster,condensed_distance_matrix):
        """
        Returns the sum of the distances of one element vs one cluster.
        """
        cdef double distance_sum
        cdef int i, N 
        cdef list all_elements
        
        all_elements =  cluster.all_elements
        N = len(all_elements)
        distance_sum = 0.0
        for i in range(N):
            distance_sum += condensed_distance_matrix[all_elements[i],element]
        return distance_sum
    
    @cython.boundscheck(False)
    cpdef double __get_average_distance_with_my_cluster(self,int element, cluster, condensed_distance_matrix):
        """
        Calculates the average distance of one element to all the elements of his cluster.
        """
        cdef double distance_sum_1 = 0
        cdef int cluster_size
        
        cluster_size = cluster.get_size()
        if cluster_size == 1:
            return 1. 
        else:
            distance_sum = self.__sum_cluster_distances(element, cluster, condensed_distance_matrix)
            return distance_sum_1 / (cluster_size-1)
    
    @cython.boundscheck(False)
    cpdef double __get_average_distance_with_cluster(self,int element,cluster,condensed_distance_matrix):
        """
        Calculates the average distance of one element to all the elements of his cluster.
        """
        cdef double distance_sum_2 = 0
        cdef int cluster_size_1
        
        cluster_size_1 = cluster.get_size()
        distance_sum = self.__sum_cluster_distances(element,cluster,condensed_distance_matrix)
        return distance_sum_2 / cluster_size_1
