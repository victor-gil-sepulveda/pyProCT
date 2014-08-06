"""
Created on 10/02/2014

@author: victor
"""
import os.path
import pyproct.protocol.saveTools as saveTools

class RepresentativesPostAction(object):

    KEYWORD = "representatives"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):
        if RepresentativesPostAction.KEYWORD in postprocessing_parameters:
            save_representatives(clustering, postprocessing_parameters[RepresentativesPostAction.KEYWORD],
                                 trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles)


def save_representatives(clustering, my_params, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):

    #Parameters
    keep_remarks = my_params.get_value("keep_remarks", default_value = False)
    keep_frame_number = my_params.get_value("keep_frame_number", default_value = False)

    # The real job
    medoids = clustering.get_medoids(matrixHandler.distance_matrix)
    # Set prototypes and ids (medoids are ordered) Once refactored this won't be necessary
    for i in range(len(clustering.clusters)):
        clustering.clusters[i].prototype = medoids[i]

    representatives_path = saveTools.save_representatives(medoids,
                                                          "representatives",
                                                          workspaceHandler,
                                                          trajectoryHandler,
                                                          do_merged_files_have_correlative_models=True,
                                                          write_frame_number_instead_of_correlative_model_number = keep_frame_number,
                                                          keep_remarks = keep_remarks )

    generatedFiles.append({"description":"Cluster representatives",
                                "path":os.path.abspath(representatives_path),
                                "type":"pdb"})

