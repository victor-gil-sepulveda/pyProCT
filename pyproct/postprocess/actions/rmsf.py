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

    def run(self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):
        try:
            displacements_path, CA_mean_square_displacements = calculate_RMSF(clustering,
                                                                                                 trajectoryHandler,
                                                                                                 workspaceHandler,
                                                                                                 matrixHandler)

            generatedFiles.append({
                                        "description":"Alpha Carbon mean square displacements",
                                        "path":os.path.abspath(displacements_path),
                                        "type":"text"
            })

            open(displacements_path,"w").write(json.dumps(CA_mean_square_displacements,
                                                  sort_keys=False,
                                                  indent=4,
                                                  separators=(',', ': ')))
        except Exception:
            print "[ERROR][Driver::postprocess] Impossible to calculate rmsf file."


def calculate_RMSF(best_clustering, trajectoryHandler, workspaceHandler, matrixHandler):
    global_cluster = Cluster(None, best_clustering.get_all_clustered_elements())
    global_cluster.prototype = global_cluster.calculate_medoid(matrixHandler.distance_matrix)
    ca_pdb_coordsets =numpy.copy(trajectoryHandler.getMergedStructure().select("name CA").getCoordsets())
    calculator = RMSDCalculator(calculatorType = "QTRFIT_SERIAL_CALCULATOR",
                                    fittingCoordsets = ca_pdb_coordsets)
    calculator.iterativeSuperposition()
    CA_mean_square_displacements= {
                                   "global":list(calc_rmsf_of_cluster(ca_pdb_coordsets,global_cluster))
                                   }
    clusters = best_clustering.clusters
    for cluster in range(len(clusters)):
        CA_mean_square_displacements[cluster.id] = superpose_and_calc_rmsf(ca_pdb_coordsets, cluster)
    displacements_path = os.path.join(workspaceHandler["results"], "CA_displacements.json")

    return displacements_path, CA_mean_square_displacements

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
