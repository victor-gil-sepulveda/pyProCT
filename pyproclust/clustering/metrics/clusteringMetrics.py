'''
Created on 30/03/2012

@author: victor
'''
import numpy
import random

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
    
    def evaluate(self,clustering,subsampling_percent,condensed_matrix):
        mean_min_dists = 0
        total_elements_in_mean_min_dist = 0
        clustering_size = len(clustering.clusters)
        for i in range(clustering_size-1):
            cluster1 = clustering.clusters[i]
            for j in range(i+1,clustering_size):
                cluster2 = clustering.clusters[j]
                (imd,jmd) = self.subsampled_mean_min_dist(cluster1, cluster2, subsampling_percent,condensed_matrix)
#                mean_min_dists.append(imd)
#                mean_min_dists.append(jmd)
                total_elements_in_mean_min_dist += 2
                mean_min_dists = mean_min_dists + imd + jmd
#        return mean_function(mean_min_dists) 
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
        
#    def calculate_mean(self,cluster1,cluster2,distance_matrix):
#        """
#        Calculates the mean of all distances between the elements of 2 clusters.
#        """
#        distances = []
#        for ei,ej in itertools.product(cluster1.all_elements,cluster2.all_elements):
#                distances.append( distance_matrix[ei,ej])
#        return mean_function(distances)

class CohesionCalculator(object):
    def __init__(self):
        pass
    
    def evaluate(self,cluster,condensed_distance_matrix):
        """
        Returns the cohesion value of a cluster. The condensed matrix given as 
        parameter stores the distances of the elements of the dataset used to extract
        the cluster.
        """
        if cluster.prototype == None:
            return self.__noproto_eval(cluster,condensed_distance_matrix)
        else:
            return self.__proto_eval(cluster,condensed_distance_matrix)
    
    def __noproto_eval(self,cluster,condensed_distance_matrix):
        """
        Evaluation of the cohesion value for a cluster without prototype.
        
        The definition of cohesion would be weight*2*cohesion in the case we follow
        the exact formula in the book []. As we are going to do comparisons, the x2 global
        multiplication doesn't affect.
        
        Cohesion of a cluster of 1 element should be infinite instead of 0...
        """
        size = cluster.get_size()
        weight = 1. / size
        cohesion = 0.
        for i in range(size-1):
            for j in range(i+1,size):
                cohesion = cohesion + condensed_distance_matrix[cluster[i],cluster[j]]
        return weight*cohesion
    
    def __proto_eval(self,cluster,condensed_distance_matrix):
        """
        Evaluation of the cohesion value for a cluster without protorype.
        """
        print "CohesionCalculator,__proto_eval Not implemented"
        exit(-1)
        
class CohesionAndSeparationCalculator(object):
    
    def __init__(self):
        pass
    
    def evaluate(self,cluster,clusterization,cluster_cohesion,condensed_distance_matrix):
        """
        Returns the cohesion plus separation value of a cluster. The condensed matrix 
        given as parameter stores the distances of the elements of the dataset used to 
        extract the cluster.
        """
        if cluster.prototype == None:
            return self.__noproto_eval(cluster,clusterization,cluster_cohesion,condensed_distance_matrix)
        else:
            return self.__proto_eval(cluster,clusterization,cluster_cohesion,condensed_distance_matrix)
        
    def __noproto_eval(self,cluster,clusterization,cluster_cohesion,condensed_distance_matrix):
        """
        Does the actual calculation for clusters without prototype.
        """
        if cluster_cohesion > 0:
            weight = 1./cluster_cohesion
            sep_and_cohe = 0.0
            ## I'm inside?
            where_am_i = clusterization.cluster_index(cluster)
            
            for i in range(len(clusterization.clusters)):    
                if i != where_am_i :
                    c_j = clusterization.clusters[i]
                    sep_and_cohe = sep_and_cohe + self.__clusters_mixed_cohesion_wo_prot(cluster,c_j,condensed_distance_matrix)
            return weight*sep_and_cohe
        else:
            return numpy.finfo(numpy.float32).max
    
    def __proto_eval(self,cluster,clusterization,cluster_cohesion,condensed_distance_matrix):
        print "CohesionAndSeparationCalculator,__proto_eval Not implemented"
        exit(-1)
        
    def __clusters_mixed_cohesion_wo_prot(self,cluster_1,cluster_2,condensed_distance_matrix):
        """
        Calculates the 'cohesion' of one cluster vs other.
        Precondition: Clusters don't have shared elements.
        """
        mixed_cohesion = 0
        for c_i in cluster_1.all_elements:
            for c_j in cluster_2.all_elements:
                mixed_cohesion = mixed_cohesion + condensed_distance_matrix[c_i,c_j]
        return mixed_cohesion

