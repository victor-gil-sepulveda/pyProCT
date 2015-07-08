'''
Created on 20/02/2014

@author: victor
'''
from pyproct.postprocess.actions.confSpaceComparison.comparator import Separator,\
    Analyzer
import os
import json

class ConfSpaceComparisonPostAction(object):
    KEYWORD = "conformational_space_comparison"

    def __init__(self):
        pass

    def run( self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, 
            matrixHandler, generatedFiles):
        
        comparison = conformational_space_comparison(clustering, trajectoryHandler,  matrixHandler, 
                         None, None)
        
        file_name = postprocessing_parameters.get_value("file", default_value = "conf_space_comp") + ".json"
        
        result_file_path = os.path.join(workspaceHandler["results"], 
                                            file_name)
        open(result_file_path, "w").write(
            json.dumps(comparison, sort_keys = False, indent = 4, separators = (',', ': '))
        )
        
        generatedFiles.append({
                               "description":"Conformational Space Comparison",
                               "path":os.path.abspath(result_file_path),
                               "type":"text"
        })
        
def conformational_space_comparison(clustering, trajectoryHandler,  matrixHandler, 
                         clustering_parameters, refinement_parameters):

#     clustering = Refiner(matrixHandler,
#                          trajectoryHandler,
#                          clustering_parameters,
#                          refinement_parameters,
#                          observer).run(clustering)

    traj_ranges = {}
    current = 0
    for i, pdb_source in enumerate(trajectoryHandler.sources):
        num_confs = pdb_source.get_info("number_of_conformations")
        traj_ranges["traj_%d"%i] = (current, current + num_confs -1)
        current = current + num_confs

    decomposed_clusters = Separator.separate(clustering.clusters, traj_ranges)

    analysis = Analyzer.run(decomposed_clusters, matrixHandler.distance_matrix)

    return analysis
