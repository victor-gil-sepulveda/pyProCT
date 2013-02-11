'''
Created on 17/09/2012

@author: victor
'''
from pyRMSD.condensedMatrix import CondensedMatrix
import pyproclust.tools.scriptTools as scripts_common
from pyproclust.clustering.selection.bestClusteringSelector import BestClusteringSelector
from pyproclust.clustering.cluster import Cluster
from pyproclust.clustering.comparison.comparator import Separator
from pyproclust.protocol.clusteringExplorationFunctions import do_clustering_exploration
from pyproclust.protocol.protocolImplementationFunctions import do_clustering_filtering, clustering_scoring
from pyRMSD.matrixHandler import MatrixHandler
import numpy

def pick_all_elements_from_clusters(clusters):
    """
    Puts all the elements of a list of clusters into a list.
    """
    all_elements = []
    for c in clusters:
        all_elements.extend(c.all_elements)
    return all_elements

def redefine_clusters_with_map(clusters,elements_map):
    """
    It renames the elements of a list of clusters using a map array.
    """
    new_clusters = []
    for cluster in clusters:
        new_cluster_elems = []
        for element in cluster.all_elements:
            new_cluster_elems.append(elements_map[element])
            
        new_clusters.append(Cluster(None, new_cluster_elems))
    return new_clusters

def recreate_matrix(matrix,number_of_elements, all_elements_map):
    """
    Picks a portion of a distance matrix and returns it as another matrix.
    """
    c_matrix = CondensedMatrix([0.0]*int(number_of_elements*(number_of_elements-1)/2.))
    for i in range(number_of_elements):
        ele_i = all_elements_map[i]
        for j in range(number_of_elements):
            ele_j = all_elements_map[j]
            c_matrix[i,j] = matrix[ele_i,ele_j]
    return c_matrix

class pureRefinementProtocol(object):
    """
    This guy tries to refine the clusters from a previous clustering.
    """
    def __init__(self,clusters,matrix,clusters_path,pdb_structure):
        # Workspace
        scripts_common.create_directory(clusters_path)
        # Pick all the elements
        self.all_elements = pick_all_elements_from_clusters(clusters)
        self.number_of_elements = len(self.all_elements)
        #Create a map
        self.all_elements_map = []
        for a in self.all_elements:
            self.all_elements_map.append(a)
        #Now recreate its matrix
        self.condensed_matrix = recreate_matrix(matrix, self.number_of_elements, self.all_elements_map)
        handler = MatrixHandler()   
        handler.distance_matrix = self.condensed_matrix
        handler.saveMatrix("matrix_backup")
        self.max_matrix_value = numpy.max(self.condensed_matrix.get_data())
        self.mean_matrix_value = numpy.mean(self.condensed_matrix.get_data())
        self.clusters_path = clusters_path
        self.pdb_structure = pdb_structure
         
    def run(self,protocol_params):
        new_clusters = self.refinement_protocol_workflow(protocol_params)
        if new_clusters != [] :
            return new_clusters
        else:
            return None
    
    def refinement_protocol_workflow(self,protocol_params):
        #OK
        do_clustering_exploration(protocol_params, get_algorithm_scheduler(protocol_params),\
                                  self.condensed_matrix, self.max_matrix_value, self.mean_matrix_value,
                                  self.clusters_path)
        #OK
        non_filtered_clusterings = scripts_common.load_binary_clusters(self.clusters_path)
        
        #OK
        filtered_clusterings = do_clustering_filtering(non_filtered_clusterings,protocol_params,\
                                                       non_filtered_clusterings,self.number_of_elements)
        new_clusters = []
        if len(filtered_clusterings) != 0:
            #OK
            string_results, results_pack = clustering_scoring(filtered_clusterings,protocol_params,\
                                                                   self.condensed_matrix,self.pdb_structure)
            print string_results
            #OK
            bestClusterSelector = BestClusteringSelector(protocol_params.cluster_score_value_map)
            best_score, best_clustering = bestClusterSelector.chooseBestClustering(results_pack) #@UnusedVariable
            
            # We need to recover the real numbering of the elements
            new_clusters = redefine_clusters_with_map(best_clustering.clusters,self.all_elements_map)
        
        return new_clusters
        
class mixedRefinementProtocol(pureRefinementProtocol):
    """
    This guy tries to refine the mixed clusters from a previous clustering. The difference
    with the first one is that in this case it can produce some pure clusters.
    """
    def __init__(self,clusters,matrix,clusters_path,pdb_structure):
        only_mixed = []
        for pack in clusters:
            only_mixed.append(pack[0])
        super( mixedRefinementProtocol, self ).__init__(only_mixed,matrix,clusters_path,pdb_structure)
        
    def run(self,protocol_params):
        new_clusters = self.refinement_protocol_workflow(protocol_params)
        
        if new_clusters != []:
            # Separate this guys
            separator = Separator()
            pure_A , pure_B, mixed = separator.c_separate(new_clusters, protocol_params.pdb1, protocol_params.pdb2)
            
            # Extract mixed clusters
            only_mixed = []
            for pack in mixed:
                only_mixed.append(pack[0])
                
            return pure_A, pure_B, only_mixed
        else:
            return None
        