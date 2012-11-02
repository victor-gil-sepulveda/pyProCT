#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
'''
Created on 12/07/2012

@author: victor
'''

import numpy
cimport numpy
import math
from libc.stdlib cimport malloc, free
cimport cython
import random

numpy.import_array()

cdef class CythonBoundedCohesionCalculator(object):
    def __init__(self):
        pass
    
    @cython.boundscheck(False)
    cpdef double evaluate(self,clustering,condensed_distance_matrix)  except *:
        cdef double weight,cohesion,max_cohesion
        cdef int size, i, j, total_number_of_elements
        
        total_number_of_elements = clustering.total_number_of_elements

        if clustering.total_number_of_elements > 0:
            max_cohesion = numpy.sum(condensed_distance_matrix.get_data())/total_number_of_elements
            total_cohesion = 0
            for c in clustering.clusters:
                size = c.get_size()
                weight = 1. / size
                cohesion = 0.
                for i in range(size-1):
                    for j in range(i+1,size):
                        cohesion = cohesion + condensed_distance_matrix[c[i],c[j]]
                total_cohesion +=  weight*cohesion
            return total_cohesion / max_cohesion
        else:
            return 0.

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
        
        # Indetermined case
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
    cpdef double __get_average_distance_with_my_cluster(self,int element,cluster,condensed_distance_matrix):
        """
        Calculates the average distance of one element to all the elements of his cluster.
        """
        cdef double distance_sum
        cdef int cluster_size
        
        cluster_size = cluster.get_size()
        if cluster_size == 1:
            return 1. 
        else:
            distance_sum = self.__sum_cluster_distances(element,cluster,condensed_distance_matrix)
            return distance_sum / (cluster_size-1)
    
    @cython.boundscheck(False)
    cpdef double __get_average_distance_with_cluster(self,int element,cluster,condensed_distance_matrix):
        """
        Calculates the average distance of one element to all the elements of his cluster.
        """
        cdef double distance_sum
        cdef int cluster_size
        
        cluster_size = cluster.get_size()
        distance_sum = self.__sum_cluster_distances(element,cluster,condensed_distance_matrix)
        return distance_sum / cluster_size

def mean_function(my_list):
    """
    Function wrapper to handle difficult cases (like empty lists :/ ).
    """
    if(len(my_list)==0):
        return 0
    else:
        return numpy.mean(my_list)   


cdef class CythonMeanMinimumDistanceCalculator(object):
    def __init__(self,seed_num = None):
        if seed_num:
            random.seed(seed_num)
    
    def evaluate(self,clustering,int subsampling_percent,condensed_matrix):
        cdef int i , j , clustering_size,total_elements_in_mean_min_dist
        cdef double mean_min_dists
        mean_min_dists = 0.
        clustering_size = len(clustering.clusters)
        total_elements_in_mean_min_dist = 0
        for i in range(clustering_size-1):
            cluster1 = clustering.clusters[i]
            for j in range(i+1,clustering_size):
                cluster2 = clustering.clusters[j]
                (imd,jmd) = self.subsampled_mean_min_dist(cluster1, cluster2, subsampling_percent,condensed_matrix)
                total_elements_in_mean_min_dist += 2
                mean_min_dists = mean_min_dists + imd + jmd
        return mean_min_dists / total_elements_in_mean_min_dist
    
    def subsampled_mean_min_dist(self,cluster1, cluster2, int subsampling_percent,condensed_matrix):
        """
        Does the calculation for two clusters. This implies that for each pair of clusters it
        gets the min_dists and those of min_dists that are smaller than the mean to get the 
        subsampled value.
        """
        min_dists, mean = self.get_mean_and_min_distances(cluster1, cluster2, condensed_matrix)
        min_dists_low_mean = self.get_distances_less_than_mean(min_dists,mean)
        
        sb1 = self.subsample(len(cluster1.all_elements), subsampling_percent, min_dists_low_mean)
        sb2 = self.subsample(len(cluster2.all_elements), subsampling_percent, min_dists_low_mean)
        
        return sb1,sb2
        
    def get_mean_and_min_distances(self,cluster1,cluster2,condensed_matrix):
        """
        Returns the minimum distances for the elements of cluster1 vs cluster2.
        """
        cdef double all_dists_accum
        min_dists = []
        all_dists_accum = 0
        for ei in cluster1.all_elements:
            distances = []
            for ej in cluster2.all_elements:
                distances.append(condensed_matrix[ei,ej])
            min_dists.append(numpy.min(distances))
            all_dists_accum += numpy.sum(distances)
        return min_dists, all_dists_accum / float(len(min_dists))
    
    cdef inline get_distances_less_than_mean(self,distances, double mean):
        """
        Returns an array with distances wich are smaller than the given mean.
        """
        a = numpy.array(distances)
        return  a[a<=mean]  
                   
    cdef double subsample(self,cluster_size,int subsampling_percent,distances):
        """
        It chooses a percent of random given distances and calculates the mean.
        """
        subsampled_elems = max(1,int(cluster_size*subsampling_percent/100.)) # minimum is 1
        random.shuffle(distances)
        return mean_function(distances[:subsampled_elems])
