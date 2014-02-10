'''
Created on 10/02/2014

@author: victor
'''
import os.path
from pyproct.tools.pdbTools import extract_frames_from_trajectory_sequentially,\
    get_number_of_frames
import pyproct.protocol.saveTools as saveTools
import tarfile




def save_all_clusters(parameters, clustering, generatedFiles, timer):
    timer.start("Save clusters")

    #Parameters
    my_params = parameters["global"]["postprocess"]["pdb_clusters"]
    keep_remarks = my_params["keep_remarks"] if "keep_remarks" in my_params else False
    keep_frame_number = my_params["keep_frame_number"] if "keep_frame_number" in my_params else False


    results_place = parameters["workspace"]["results"]
    clusterings_place = parameters["workspace"]["clusterings"]
    tmp_place = parameters["workspace"]["tmp"]

    # The real job
    input_path = os.path.join(tmp_place, "tmp_merged_trajectory.pdb")
    number_of_frames = get_number_of_frames(input_path)
    bz2_files = []
    for cluster in clustering:
        output_path = os.path.join(clusterings_place, "%.bz2"%(cluster.id))
        bz2_files.append(output_path)

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
    tar_path = os.path.join(results_place,"clusters.tar.bz2")
    tar = tarfile.open(tar_path, "w:bz2")
    for comp_file in bz2_files:
        tar.add(comp_file)
    tar.close()
    timer.stop("Save clusters")

    generatedFiles.append({"description":"Frames by cluster",
                                         "path":tar_path,
                                         "type":"compressed_pdb"})()

def save_representatives(clustering, parameters, matrixHandler, workspaceHandler,
                         trajectoryHandler, generatedFiles, timer):
    timer.start("Representatives")

    #Parameters
    my_params = parameters["global"]["postprocess"]["pdb_clusters"]
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

    generatedFiles.append({"description":"Cluster central conformations (II)",
                                "path":representatives_path,
                                "type":"pdb"})

    timer.stop("Representatives")
