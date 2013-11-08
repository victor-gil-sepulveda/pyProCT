'''
Created on 16/04/2012

@author: victor
'''
from pyproclust.clustering.cluster import gen_clusters_from_class_list
from pyproclust.clustering.clustering import Clustering


class PointClassType:
    UNCLASSIFIED = -1
    NOISE = 0

class DBSCANAlgorithm(object):
    '''
    Implementation of the DBSCAN algorithm (Ester et al. KDD-96)
    '''

    def __init__(self,condensed_matrix, **kwargs):
        self.condensed_matrix = condensed_matrix
        self.number_of_elements = self.condensed_matrix.get_number_of_rows()
        self.element_class = [PointClassType.UNCLASSIFIED]*self.number_of_elements
        
    def perform_clustering(self,kwargs):
        """
        Main loop to perform the DBSCAN algorithm. 
        """
        eps = kwargs["eps"]
        minpts = kwargs["minpts"]
        current_cluster_id = PointClassType.NOISE + 1
        for i in range(self.number_of_elements):
            current_element = i
            if self.element_class[current_element] == PointClassType.UNCLASSIFIED:
                if self.__expand_cluster(current_element,current_cluster_id,eps,minpts):
                    current_cluster_id = current_cluster_id + 1
    
        # Return the clusters once the clustering is done
        # NOISE elements form a single cluster with ID = PointClassType.NOISE 
        # and will be removed from the clustering
        clusters = gen_clusters_from_class_list(self.element_class,skip_list = [PointClassType.NOISE])
        return Clustering(clusters,details="DBSCAN (eps = "+str(eps)+" minpts = "+str(minpts)+") "+str(self.number_of_elements)+" elems")
    
    def __expand_cluster(self, current_element, current_cluster_id, eps, minpts):
        """
        Secondary loop.
        """
#         seeds = self.__eps_neighborhood(current_element,eps)
        seeds = self.condensed_matrix.element_neighbors_within_radius(current_element,eps)
        if len(seeds) < minpts:
            self.element_class[current_element] = PointClassType.NOISE
            return False
        else:
            self.element_class[current_element] = current_cluster_id
            for s in seeds:
                self.element_class[s] = current_cluster_id
            self.__seed_expansion(current_cluster_id, eps, minpts, seeds)
            return True        
    
    def __seed_expansion(self, current_cluster_id, eps, minpts, seeds):
        """
        Performs iteratively a graph generation algorithm. Two nodes are connected
        if the node is a core node and the other is in its eps-neighborhood.
        """
#         while len(seeds) > 0:
#             current = seeds.pop()
#             extra_seeds = self.__eps_neighborhood(current, eps)
#             if len(extra_seeds) >= minpts:
#                 for s in extra_seeds:
#                     if self.element_class[s] in [PointClassType.UNCLASSIFIED, PointClassType.NOISE]:
#                         if self.element_class[s] == PointClassType.UNCLASSIFIED:
#                             seeds.append(s)
#                         self.element_class[s] = current_cluster_id
        seeds = list(seeds)
        while len(seeds) > 0:
            current = seeds.pop()
            extra_seeds = self.condensed_matrix.element_neighbors_within_radius(current, eps)
            if len(extra_seeds) >= minpts:
                for s in extra_seeds:
                    if self.element_class[s] in [PointClassType.UNCLASSIFIED, PointClassType.NOISE]:
                        if self.element_class[s] == PointClassType.UNCLASSIFIED:
                            seeds.append(s)
                        self.element_class[s] = current_cluster_id
            
    # __eps_neighborhood will not return the central element!!
    def __eps_neighborhood(self,node,eps):
        """
        Populates the neighbor list for a given node (without the node itself).
        """
        neighbours = []
        for i in range(self.number_of_elements):
            if node != i and self.condensed_matrix[node,i] <= eps:
                neighbours.append(i)
        return neighbours
        