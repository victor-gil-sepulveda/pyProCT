'''
Created on 07/08/2012

@author: victor
'''

class ClusteringFilter(object):
    """
    This class implements functions to get a list of clusterings and return a list
    of filtered clusterings (hopefully smaller).
    """
    
    def __init__(self, evaluation_parameters, matrix_handler):
        """
        Constructor.
        
        @param evaluation_parameters: Structure that at least defines: 
            - evaluation_parameters['minimum_clusters'] : Minimum number of clusters allowed for a clustering.
            - evaluation_parameters['maximum_clusters'] : Maximum number of clusters allowed for a clustering.
            - evaluation_parameters['minimum_cluster_size'] : Minimum number of elements per cluster (clusters
              with less than this will be treated as noise).
            - evaluation_parameters['maximum_noise'] : Max percentage of noise allowed (non clustered elements).
            
        @param distance_matrix: Is the current distance matrix.
        """
        self.evaluation_parameters = evaluation_parameters
        
        #'total_number_of_elements' is the total number of elements of the dataset (which may coincide or not
        #with the number of clustered elements).
        self.total_number_of_elements = matrix_handler.distance_matrix.row_length
    
    def filter(self, clusterings):
        """
        Eliminates the clusters whose parameters are not in the range we have defined in the script evaluation parameters. 
        
        @param clusterings: Is a clustering_info dictionary indexed by clustering_id
        
        @return: A tuple containing the clustering_info structures of the selected and not selected clusterings. These last
        with their reasons for not having been selected.
        """
        selected_clusterings =  {}
        not_selected_clusterings = {}
        
        for clustering_id in clusterings:
            reasons = self.check_clustering(clusterings[clustering_id]["clustering"])
            
            if(reasons == []):
                selected_clusterings[clustering_id] = clusterings[clustering_id]
            else:
                not_selected_clusterings[clustering_id] = clusterings[clustering_id]
                not_selected_clusterings[clustering_id]["reasons"] = reasons
        
        return selected_clusterings, not_selected_clusterings
    
    def check_num_clusters_in_range(self, clustering):
        """
        Used to see if a clustering has a number of clusters in the range defined by script's evaluation parameters.
        
        @param clustering: The clustering to be checked.
        
        @return: The reasons of not being selected (if any).
        """
        num_clusters = len(clustering.clusters)
        
        if num_clusters <  self.evaluation_parameters["minimum_clusters"]:
            return [{"reason":"TOO_FEW_CLUSTERS","data":{"minimum":self.evaluation_parameters["minimum_clusters"],"current":num_clusters}}]
        
        if num_clusters >  self.evaluation_parameters["maximum_clusters"]:
            return [{"reason":"TOO_MUCH_CLUSTERS","data":{"maximum":self.evaluation_parameters["maximum_clusters"],"current":num_clusters}}]
       
        return []
    
    def check_noise_level(self, clustering):
        """
        Does the noise level check.
        
        @param clustering: The clustering to be checked.
        
        @return: The reasons of not being selected (if any).
        """
        noise_level = 100 - (float(clustering.total_number_of_elements) / self.total_number_of_elements) *  100
        
        if noise_level > self.evaluation_parameters['maximum_noise']:
            return  [{"reason":"TOO_MUCH_NOISE","data":{"maximum":self.evaluation_parameters["maximum_noise"],"current":noise_level}}]
        else:
            return []
        
    def check_clustering(self, clustering):
        """
        Decides if one clustering is inside the defined parameters, and returns the reasons of why not otherwise.
        
        @param clustering: The clustering to be checked.
        
        @return: The reasons of not being selected (if any).
        """
        # First we must eliminate the noise
        clustering.eliminate_noise(self.evaluation_parameters["minimum_cluster_size"])
        
        # Then, look for reasons to reject this clustering
        reasons = self.check_num_clusters_in_range(clustering)
        reasons.extend(self.check_noise_level(clustering))

        return reasons 
