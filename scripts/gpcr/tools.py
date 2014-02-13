import json
from pyproct.tools.commonTools import convert_to_utf8
import numpy
import os


PYPROCT = "/home/victor/workspaces/Python/pyProClust/pyproct/main.py"

def load_dic_in_json(filename):
    return convert_to_utf8(json.loads("".join(open(filename,"r").readlines())))

def save_dic_in_json(this, there):
    open(there, "w").write(json.dumps(this, sort_keys=False, indent=4, separators=(',', ': ')))

def create_dir(_dir):
    if not os.path.exists(_dir):
        os.makedirs(_dir)

def use_pyproct(directory, script):
    create_dir(directory)
    os.system("python %s %s "%(PYPROCT, script))

def get_best_clustering(clustering_directory):
    results = load_dic_in_json(os.path.join(clustering_directory,"results","results.json"))
    best_clustering_id = results['best_clustering']
    return results['selected'][best_clustering_id]

def score_cluster(cluster, matrix, metrics):
    all_rmsd_distances = [metrics[element][0] for element in cluster.all_elements]
    all_binding_energies = [metrics[element][1] for element in cluster.all_elements]
    return numpy.mean(all_rmsd_distances),numpy.mean(all_binding_energies)

def find_most_negative_be_half(scores):
    # Retrieve all scores and do the mean
    be_scores_mean = numpy.mean([x[0][1] for x in scores])
    filtered_scores = []
    for score in scores:
        rmsd,be = score[0]
        if be <= be_scores_mean:
            filtered_scores.append(score)
    return filtered_scores

def find_most_negative_rmsd_half(scores):
    # Retrieve all scores and do the mean
    rmsd_scores_mean = numpy.mean([x[0][0] for x in scores])
    filtered_scores = []
    for score in scores:
        rmsd,be = score[0]
        if rmsd <= rmsd_scores_mean:
            filtered_scores.append(score)
    return filtered_scores

def find_most_negative_cluster(scores):
    filtered = scores
    while len(filtered)>1:
        filtered = find_most_negative_be_half(filtered)
        filtered = find_most_negative_rmsd_half(filtered)
    return filtered[0]

def find_5_clusters_with_less_energy(scores):
    # Mean energy will be the second element of each score
    energies = [(en,cluster) for ((dist, en),cluster) in scores]  # @UnusedVariable
    energies.sort()
    clusters = []
    for i in range(5):
        clusters.append(energies[i])
    return clusters

def normalize_metrics(metrics):
    normalized_metrics = []
    metricsT = metrics.T
    for metric_array in metricsT:
        m_max = numpy.max(metric_array)
        m_min = numpy.min(metric_array)
        normalized_metrics.append( (metric_array - m_min) / (m_max - m_min))
    return numpy.array(normalized_metrics).T

