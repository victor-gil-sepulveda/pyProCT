'''
Created on 20/02/2014

@author: victor
'''
from pyproct.protocol.refinement.Refiner import Refiner
from pyproct.clustering.comparison.comparator import Separator, Analyzer


def conformational_space_comparison(clustering, matrixHandler, trajectoryHandler,
                         clustering_parameters, refinement_parameters,  observer):

#     clustering = Refiner(matrixHandler,
#                          trajectoryHandler,
#                          clustering_parameters,
#                          refinement_parameters,
#                          observer).run(clustering)

    # TODO: testing
    traj_ranges = {}
    current = 0
    for i, pdb in enumerate(trajectoryHandler.pdbs):
        traj_ranges["traj_%d"%i] = (current, current + pdb["conformations"] -1)
        current = current + pdb["conformations"]

    decomposed_clusters = Separator.separate(clustering.clusters, traj_ranges)

    analysis = Analyzer.run(decomposed_clusters, matrixHandler.distance_matrix)

    return analysis
