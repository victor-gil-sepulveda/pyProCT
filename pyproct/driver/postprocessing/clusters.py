'''
Created on 10/02/2014

@author: victor
'''
import os.path
from pyproct.tools.pdbTools import extract_frames_from_trajectory_sequentially,\
    get_number_of_frames
import pyproct.protocol.saveTools as saveTools
import tarfile
import bz2
from pyproct.protocol.saveTools import merge_pdbs

def save_all_clusters(my_params, pdb_params, workspaceHandler, clustering, generatedFiles, timer):
    timer.start("Save clusters")

    #Parameters
    keep_remarks = my_params["keep_remarks"] if "keep_remarks" in my_params else False
    keep_frame_number = my_params["keep_frame_number"] if "keep_frame_number" in my_params else False

    # Places
    results_place = workspaceHandler["results"]
    clusters_place = workspaceHandler["clusters"]
    tmp_place = workspaceHandler["tmp"]

    # The real job
    input_path = os.path.join(tmp_place, "tmp_merged_trajectory.pdb")
    merge_pdbs(pdb_params, input_path)

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
    timer.stop("Save clusters")

    generatedFiles.append({"description":"Clusters",
                                         "path":tar_path,
                                         "type":"compressed_pdb"})

def save_representatives(clustering, my_params, matrixHandler, workspaceHandler,
                         trajectoryHandler, generatedFiles, timer):
    timer.start("Representatives")

    #Parameters
    keep_remarks = my_params["keep_remarks"] if "keep_remarks" in my_params else False
    keep_frame_number = my_params["keep_frame_number"] if "keep_frame_number" in my_params else False

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
                                "path":representatives_path,
                                "type":"pdb"})

    timer.stop("Representatives")
