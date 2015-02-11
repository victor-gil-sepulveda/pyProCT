"""
Created on 16/07/2014

@author: victor
"""
import json
import os
from pyRMSD.RMSDCalculator import RMSDCalculator
from pyproct.clustering.cluster import Cluster
import numpy

class RmsfPostAction(object):
    KEYWORD = "rmsf"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, data_handler, workspaceHandler, matrixHandler, generatedFiles):
        
        rmsf_per_cluster = calculate_RMSF(clustering, data_handler)

        rmsf_file_path = os.path.join(workspaceHandler["results"], "CA_displacements.json")
        open(rmsf_file_path,"w").write(json.dumps(rmsf_per_cluster,
                                              sort_keys=False,
                                              indent=4,
                                              separators=(',', ': ')))
        
        generatedFiles.append({
                                "description":"Alpha Carbon mean square displacements",
                                "path":os.path.abspath(rmsf_file_path),
                                "type":"text"
        })

def calculate_RMSF(best_clustering, data_handler):
    ca_pdb_coordsets = data_handler.get_data().getSelectionCoordinates("name CA")

    global_cluster = Cluster(None, best_clustering.get_all_clustered_elements())
    global_cluster.id = "global"
    
    clusters = best_clustering.clusters + [global_cluster]
    rmsf_per_cluster = {}
    for cluster in clusters:
        rmsf_per_cluster[cluster.id] = superpose_and_calc_rmsf(ca_pdb_coordsets, cluster)
    

    return rmsf_per_cluster

def superpose_and_calc_rmsf(ca_pdb_coordsets, cluster):
    # Pick the coordinates (ensuring that we are copying them)
    fitting_coordinates_of_this_cluster = ca_pdb_coordsets[cluster.all_elements]

    calculator = RMSDCalculator(calculatorType = "QTRFIT_SERIAL_CALCULATOR",
                                fittingCoordsets = fitting_coordinates_of_this_cluster)

    # Make an iterative superposition (to get the minimum RMSD of all with respect to a mean conformation)
    calculator.iterativeSuperposition()
    
    return list(calc_rmsf_of_cluster(fitting_coordinates_of_this_cluster, cluster))

def calc_rmsf_of_cluster(cluster_coordsets, cluster):
    if len(cluster_coordsets) != 0:
        mean_conformation = cluster_coordsets.mean(0)
        ssqf = numpy.zeros(mean_conformation.shape)
        for conf in cluster_coordsets:
            ssqf += (conf - mean_conformation) ** 2
        return (ssqf.sum(1) / cluster_coordsets.shape[0])**0.5
    else:
        print "[WARNING][calc_rmsf_of_cluster]  No CA atoms found. Aborting operation."
        return
