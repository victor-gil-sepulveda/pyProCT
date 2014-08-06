"""
Created on 09/01/2013

@author: victor
"""
import random
import numpy

def mean_function(my_list):
    """
    Function wrapper to handle difficult cases (like empty lists :/ ).
    """
    if(len(my_list)==0):
        return 0
    else:
        return numpy.mean(my_list)

class MeanMinimumDistanceCalculator(object):
    def __init__(self,seed_num = None):
        if seed_num:
            random.seed(seed_num)
    
    def evaluate(self, clustering, condensed_matrix, subsampling_percent = 20):
        mean_min_dists = 0
        total_elements_in_mean_min_dist = 0
        clustering_size = len(clustering.clusters)
        for i in range(clustering_size-1):
            cluster1 = clustering.clusters[i]
            for j in range(i+1,clustering_size):
                cluster2 = clustering.clusters[j]
                (imd,jmd) = self.subsampled_mean_min_dist(cluster1, cluster2, subsampling_percent,condensed_matrix)
                total_elements_in_mean_min_dist += 2
                mean_min_dists = mean_min_dists + imd + jmd
        return mean_min_dists / total_elements_in_mean_min_dist
    
    def subsampled_mean_min_dist(self,cluster1, cluster2, subsampling_percent,condensed_matrix):
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
        min_dists = []
#        all_dists = []
        all_dists_accum = 0
        for ei in cluster1.all_elements:
            distances = []
            for ej in cluster2.all_elements:
                distances.append(condensed_matrix[ei,ej])
            min_dists.append(numpy.min(distances))
#            all_dists.extend(distances)
            all_dists_accum += numpy.sum(distances)
        return min_dists, all_dists_accum / float(len(min_dists))
    
    def get_distances_less_than_mean(self,distances, mean):
        """
        Returns an array with distances wich are smaller than the given mean.
        """
        a = numpy.array(distances)
        return  a[a<=mean]  
                   
    def subsample(self,cluster_size,subsampling_percent,distances):
        """
        It chooses a percent of random given distances and calculates the mean.
        """
        subsampled_elems = max(1,int(cluster_size*subsampling_percent/100.)) # minimum is 1
        random.shuffle(distances)
        return mean_function(distances[:subsampled_elems])
