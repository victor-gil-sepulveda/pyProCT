'''
Created on 30/06/2014

@author: victor
'''
from pyproct.clustering.metrics.common import get_intra_cluster_distances
import os.path


def calculate_per_cluster_stats(best_clustering, matrix, parameters, results_folder):
    """
    CSV file
    """
    file_name = parameters.get_value("file", default_value = "per_cluster_stats") + ".csv"
    stats_file = open(os.path.join(results_folder,file_name),"w")
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



