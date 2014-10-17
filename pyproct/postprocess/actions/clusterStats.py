"""
Created on 30/06/2014

@author: victor
"""
import os.path
from pyproct.clustering.evaluation.metrics.common import get_intra_cluster_distances,\
    update_medoids, get_distances_of_elements_to
from pyproct.tools.exceptions import SingularClusterException


class ClusterStatsPostAction(object):
    KEYWORD = "cluster_stats"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, data_handler, workspaceHandler, matrixHandler, generatedFiles):
        stats_file_path = calculate_per_cluster_stats(clustering,
                                                      matrixHandler.distance_matrix,
                                                      postprocessing_parameters,
                                                      workspaceHandler["results"])
        generatedFiles.append({
                                    "description":"Stats for all clusterings (diameter and distances from center)",
                                    "path":os.path.abspath(stats_file_path),
                                    "type":"text"
        })

def calculate_per_cluster_stats(best_clustering, matrix, parameters, results_folder):
    """
    CSV file
    """
    file_name = parameters.get_value("file", default_value = "per_cluster_stats") + ".csv"
    stats_file_path = os.path.join(results_folder,file_name)
    stats_file = open(stats_file_path,"w")
    header_line =","
    for i in range(len(best_clustering.clusters)):
        cluster = best_clustering.clusters[i]
        header_line+="%s,"%cluster.id
    header_line = header_line[:-1] +"\n"

    stats_file.write(header_line)

    # TODO: Once clusterings and clusters become inmutable its medoids will be always updated,
    # then this kind of operations will be unnecessary 
    update_medoids(best_clustering, matrix)
    #----------------------------------------
    
    for i in range(len(best_clustering.clusters)):
        cluster_i = best_clustering.clusters[i]
        
        try:
            intra_distances = get_intra_cluster_distances(cluster_i, matrix)
            diameter = max(intra_distances) 
            distances_from_proto = get_distances_of_elements_to(cluster_i.prototype, 
                                                                cluster_i.all_elements, 
                                                                matrix)
            radius = max(distances_from_proto)
        except SingularClusterException:
            diameter = 0
            radius = 0
        finally:
            line = "%s(d: %.2f r: %.2f),"%(cluster_i.id, diameter, radius)

        for j in range(0, i+1):
            line += ","

        for j in range(i+1, len(best_clustering.clusters)):
            cluster_j = best_clustering.clusters[j]
            line+="%.2f,"%matrix[ cluster_i.prototype, cluster_j.prototype]

        line = line[:-1] + "\n"
        stats_file.write(line)
    stats_file.close()
    return stats_file_path



