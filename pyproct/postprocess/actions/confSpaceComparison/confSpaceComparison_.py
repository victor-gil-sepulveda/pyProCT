'''
Created on 20/02/2014

@author: victor
'''
from pyproct.postprocess.actions.confSpaceComparison.comparator import Separator,\
    Analyzer


class ConfSpaceComparisonPostAction(object):
    KEYWORD = "conf_space_comparison"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):
        pass


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
