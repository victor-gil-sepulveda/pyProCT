"""
Created on 16/07/2014

@author: victor
"""

import prody
import os.path
import tarfile
from pyproct.tools.pdbTools import extract_frames_from_trajectory_sequentially,\
    get_number_of_frames

class SaveAllClustersPostAction(object):
    KEYWORD = "representatives"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):
        save_all_clusters(clustering, postprocessing_parameters[SaveAllClustersPostAction.KEYWORD],
                          workspaceHandler, trajectoryHandler, generatedFiles)

def save_all_clusters(clustering, my_params, workspaceHandler, trajectoryHandler,  generatedFiles):

    #Parameters
    keep_remarks = my_params.get_value("keep_remarks", default_value = False)
    keep_frame_number = my_params.get_value("keep_frame_number", default_value = False)

    # Places
    results_place = workspaceHandler["results"]
    clusters_place = workspaceHandler["clusters"]
    tmp_place = workspaceHandler["tmp"]


    merged_pdb = trajectoryHandler.getMergedStructure()
    input_path = os.path.join(tmp_place, "tmp_merged_trajectory.pdb")
    prody.writePDB(input_path, merged_pdb)

    number_of_frames = get_number_of_frames(input_path)
    cluster_files = []
    for cluster in clustering.clusters:
        output_path = os.path.join(clusters_place, "%s.pdb"%(cluster.id))
        cluster_files.append(output_path)
        file_handler_in = open(input_path,"r")
        file_handler_out = open(output_path,"w")
        extract_frames_from_trajectory_sequentially(file_handler_in,
                                                    number_of_frames,
                                                    file_handler_out,
                                                    cluster.all_elements,
                                                    keep_header = keep_remarks,
                                                    write_frame_number_instead_of_correlative_model_number=keep_frame_number)
        file_handler_in.close()
        file_handler_out.close()

    # Add all bz2 files to a tar file
    tar_path = os.path.join(results_place,"clusters.tar.gz")
    tar = tarfile.open(tar_path, "w:gz")
    for comp_file in cluster_files:
        tar.add(comp_file, os.path.basename(comp_file))
    tar.close()

    generatedFiles.append({"description":"Clusters",
                                         "path":os.path.abspath(tar_path),
                                         "type":"compressed_pdb"})
