'''
Created on 07/08/2012

@author: victor
'''

class ClusteringFilter(object):
    """
    This class implements functions to get a list of clusterings and return a list
    of filtered clusterings (hopefully smaller).
    """
    def __init__(self, params):
        '''
        'params' has to define: 
            params.min_clusters - Min number of clusters allowed for a clustering
            params.max_clusters - Max number of clusters allowed for a clustering
            params.min_cluster_size - Min number of elements per cluster (clusters
                                    with less than this will be treated as noise)
            params.max_noise - Max percentage of noise allowed (non clustered elements)
        '''
        self.params = params
    
    def doClusteringFiltering(self,clusterings,total_num_elements):
        """
        Eliminates the clusters whose parameters are not in the range we have
        defined in this parameters instance. total_num_elements is the number
        of elements of the dataset.
        """
        selected_clusterings =  []
        not_selected_clusterings = []
        
        for c in clusterings:
            reasons = self.__clusteringIsAllowed(c,total_num_elements)
            if(reasons==""):
                selected_clusterings.append(c)
            else:
                not_selected_clusterings.append((c,reasons))
        return selected_clusterings, not_selected_clusterings
    
    def __numClustersInRange(self,num_clusters):
        """
        Used to see if a clustering has a number of clusters in range.
        """
        return num_clusters>=self.params.min_clusters and num_clusters<=self.params.max_clusters
    
    def __clusteringIsAllowed(self,clustering,total_num_elements):
        """
        Decides if one clustering is inside the defined params, and
        returns a reason of why not otherwise.
        'clustering' is the clustering to study and 'total_num_elements' is 
        the total number of elements of the dataset (which may coincide or not
        with the number of clustered elements).
        """
        reason = ""
        # First we must eliminate the noise
        clustering.eliminate_noise(self.params.min_cluster_size)
        
        num_clusters = len(clustering.clusters)
        num_elements = clustering.total_number_of_elements
        noise_level = 100 - (float(num_elements) / total_num_elements) *  100
        
        if not self.__numClustersInRange(num_clusters):
            reason += "Out num. clusters range:[%d, %d] with %d clusters\n"%(self.params.min_clusters,self.params.max_clusters,num_clusters)
        
        if noise_level >= self.params.max_noise:
            reason += "Too much noise (max = %.3f%%) with %.3f%% of noise\n"%(self.params.max_noise, noise_level)
        if reason != "":
            print reason
        return reason 
