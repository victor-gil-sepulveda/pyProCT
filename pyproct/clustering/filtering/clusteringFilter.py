"""
Created on 07/08/2012

@author: victor
"""

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
    
    def add_filtered_clustering(self, selected, not_selected, clusterings, clustering_id, reasons):
        """
        Adds a clustering to selected or not_selected arrays depending on if the rejection reasons
        array is empty or not.
        @param selected: The 'clustering info' dictionary with accepted clusterings.
        @param not_selected: The 'clustering info' dictionary with rejected clusterings.
        @param clusterings: The 'clustering info' dictionary.
        @param clustering_id: Is the id in a 'clustering info' dictionary of the clustering being checked.
        @param reasons: Rejection reasons array.
        """
        if(reasons == []):
            selected[clustering_id] = clusterings[clustering_id]
        else:
            not_selected[clustering_id] = clusterings[clustering_id]
            not_selected[clustering_id]["reasons"] = reasons
    
    def filter(self, clusterings):
        """
        Eliminates the clusters whose parameters are not in the range we have defined in the script evaluation parameters. 
        
        @param clusterings: Is a clustering_info dictionary indexed by clustering_id
        
        @return: A tuple containing the clustering_info structures of the selected and not selected clusterings. These last
        with their reasons for not having been selected.
        """
        selected_clusterings = {}
        not_selected_clusterings = {}
        
        # Filterings based in clustering analysis 
        for clustering_id in clusterings:
            reasons = self.check_clustering(clusterings[clustering_id]["clustering"])
            self.add_filtered_clustering(selected_clusterings,
                                         not_selected_clusterings,
                                         clusterings,
                                         clustering_id,
                                         reasons)
        
        return self.filter_repeated(selected_clusterings, not_selected_clusterings)
    
    def filter_repeated(self, selected_clusterings, not_selected_clusterings):
        """
        Checks if the clustering is already inside 'all_clusterings', returning a rejection reason in this case.
        
        @param clustering_id: The clustering id to be tested.
        @param all_clusterings: The array containing all the 'clustering info' dictionaries.
        
        @return: The new selected and not selected 'clustering info' arrays
        """
        new_selected = {}
        selected_ids = selected_clusterings.keys()
        
        for i in range(len(selected_ids)):
            clustering_id = selected_ids[i]
            clustering = selected_clusterings[clustering_id]["clustering"]
            reasons = []
            
            for j in range(i+1,len(selected_ids)):
                other_clustering_id = selected_ids[j]
                other_clustering = selected_clusterings[other_clustering_id]["clustering"]
                if clustering == other_clustering:
                    reasons = [{"reason":"EQUAL_TO_OTHER_CLUSTERING","data":{"id":other_clustering_id}}]
                    break
                
            self.add_filtered_clustering(new_selected,
                                         not_selected_clusterings,
                                         selected_clusterings,
                                         clustering_id,
                                         reasons)
             
        return new_selected, not_selected_clusterings
    
    def check_is_not_repeated(self, clustering, all_clusterings):
        """
        Checks that the clustering is not repeated again in the clustering list. Clusterings are 
        considered equal if their clusters are equal (even if they were generated using different
        algorithms and / or parameters).
        
        @param clustering: The clustering to be checked.
        
        @param all_clusterings: The array containing all the 'clustering info' dictionaries.
        
        @return: A dictionary with the reasons to discard the clustering if any.
        """
        for other_clustering in all_clusterings:
            if clustering == all_clusterings[other_clustering]:
                return [{"reason":"EQUAL_TO_OTHER_CLUSTERING","data":{"id":other_clustering}}]
        return []
    
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
