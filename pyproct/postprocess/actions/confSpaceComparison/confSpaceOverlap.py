'''
Created on Mar 4, 2016

@author: victor
'''
import os
import json
from pyproct.postprocess.actions.confSpaceComparison.comparator import Separator
import math
from scipy.stats import entropy
from numpy.linalg import norm
import numpy

import matplotlib.pyplot as plt
import seaborn as sns
 
sns.set_style("whitegrid")

class ConfSpaceOverlapPostAction(object):
    KEYWORD = "conformational_space_overlap"

    def __init__(self):
        pass

    def run( self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, 
            matrixHandler, generatedFiles):
        
        comparison = conformational_space_overlap(clustering, trajectoryHandler,  matrixHandler)
        
        file_name = postprocessing_parameters.get_value("file", default_value = "conf_space_overlap") + ".json"
        
        result_file_path = os.path.join(workspaceHandler["results"], 
                                            file_name)
        open(result_file_path, "w").write(
            json.dumps(comparison, sort_keys = False, indent = 4, separators = (',', ': '))
        )
        
        generatedFiles.append({
                               "description":"Conformational Space Overlap",
                               "path":os.path.abspath(result_file_path),
                               "type":"text"
        })
        
def conformational_space_overlap(clustering, trajectoryHandler,  matrixHandler):
    current = 0
    traj_ranges = {}
    traj_to_file = {}
    total_populations = {}
    for i, pdb_source in enumerate(trajectoryHandler.sources):
        num_confs = pdb_source.get_info("number_of_conformations")
        traj_id = "traj_%d"%i
        traj_ranges[traj_id] = (current, current + num_confs -1)
        traj_to_file[traj_id] = pdb_source.get_path()
        total_populations[traj_id] = num_confs
        current = current + num_confs

    decomposed_clusters = Separator.decompose(clustering.clusters, traj_ranges)
    
    # Get population percents for each cluster
    all_traj_ids = total_populations.keys()
    relative_populations = []
    for cluster_id in decomposed_clusters:
        dc = decomposed_clusters[cluster_id]
        relative_population =  []
        for traj_id in all_traj_ids:
            if traj_id in dc:
                relative_population.append(len(dc[traj_id]) / float(total_populations[traj_id]))
            else:
                relative_population.append(0.)
        relative_population.append(cluster_id)
        relative_populations.append(tuple(relative_population))
    
    # Sort by first traj (to 'prettify' it a bit)
    relative_populations.sort()
    cluster_ids = [rp[-1] for rp in relative_populations]
    relative_populations = numpy.array([rp[0:len(all_traj_ids)] for rp in relative_populations])
    rel_pop_per_id = {}
    sm_rel_pop_per_id = {}
    for i in range(len(all_traj_ids)):
#         print all_traj_ids[i]
#         print relative_populations.T[i]
        rel_pop_per_id[all_traj_ids[i]] = list(relative_populations.T[i])
        sm_rel_pop_per_id[all_traj_ids[i]] = smoothed(relative_populations.T[i])
        plt.plot(relative_populations.T[i],label = traj_to_file[all_traj_ids[i]])
    plt.legend()
    plt.show()
    # Calculate JSDs
    jsds = {}
    for traj_a in all_traj_ids:
        jsds[traj_a] = {}
        for traj_b in all_traj_ids:
            jsds[traj_a][traj_b] = JSD(sm_rel_pop_per_id[traj_a],
                                       sm_rel_pop_per_id[traj_b])
            
    # Compile results
    results = {
               "id_to_path":traj_to_file,
               "populations": rel_pop_per_id,
               "JSD": jsds,
               "cluster_ids": cluster_ids
    }
    
    return results
    
def smoothed(distribution, small_value = 1.0e-8):
    """
    Applies a smoothing process to the distribution.
    See http://mathoverflow.net/questions/72668/how-to-compute-kl-divergence-when-pmf-contains-0s
    for an explanation about the problem and the solution.
     
    @param distribution: distribution to be smoothed
    @param small_value: value to be set to those bins with 0 probability
     
    @return: The smoothed distribution.
    """
    total_number_of_samples = len(distribution)
    samples_in_distrib = numpy.count_nonzero(distribution)
    if samples_in_distrib > 0:
        pc = small_value * (total_number_of_samples - samples_in_distrib) / samples_in_distrib
        smoothed_distrib = numpy.empty(len(distribution))
        for i in range(len(distribution)):
            if distribution[i] == 0:
                smoothed_distrib[i] = small_value
            else:
                smoothed_distrib[i] = distribution[i] - pc
        return numpy.array(smoothed_distrib)
    else:
        return distribution

def JSD(P, Q):
    """
    Calculates the Jensen-Shannon divergence as a metric (sq_root)
    See: http://www.researchgate.net/publication/3084774_A_new_metric_for_probability_distributions
    """
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return math.sqrt(0.5 * (entropy(_P, _M) + entropy(_Q, _M)))
    
    