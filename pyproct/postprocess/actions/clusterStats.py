"""
Created on 30/06/2014

@author: victor
"""
import os.path
from pyproct.clustering.evaluation.metrics.common import get_intra_cluster_distances


class ClusterStatsPostAction(object):
    KEYWORD = "cluster_stats"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):
        stats_file_path = calculate_per_cluster_stats(clustering,
                                                      matrixHandler.distance_matrix,
                                                      postprocessing_parameters[ClusterStatsPostAction.KEYWORD],
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

    for i in range(len(best_clustering.clusters)):
        cluster_i = best_clustering.clusters[i]
        intra_distances = get_intra_cluster_distances(cluster_i, matrix)
        radius = max(intra_distances) if intra_distances != [] else 0.
        line = "%s(%.2f),"%(cluster_i.id, radius)

        for j in range(0, i+1):
            line += ","

        for j in range(i+1, len(best_clustering.clusters)):
            cluster_j = best_clustering.clusters[j]
            line+="%.2f,"%matrix[ cluster_i.prototype, cluster_j.prototype]

        line = line[:-1] + "\n"
        stats_file.write(line)
    stats_file.close()
    return stats_file_path



