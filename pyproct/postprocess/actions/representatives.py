"""
Created on 10/02/2014

@author: victor
"""
import os.path
import prody
import pyproct.tools.pdbTools as pdb_tools

class RepresentativesPostAction(object):

    KEYWORD = "representatives"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):
        save_representatives_wrapper(clustering, postprocessing_parameters[RepresentativesPostAction.KEYWORD],
                             trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles)


def save_representatives_wrapper(clustering, my_params, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):

    #Parameters
    keep_remarks = my_params.get_value("keep_remarks", default_value = False)
    keep_frame_number = my_params.get_value("keep_frame_number", default_value = False)

    # The real job
    medoids = clustering.get_medoids(matrixHandler.distance_matrix)
    # Set prototypes and ids (medoids are ordered) Once refactored this won't be necessary
    for i in range(len(clustering.clusters)):
        clustering.clusters[i].prototype = medoids[i]

    representatives_path = save_representatives(medoids,
                                                          "representatives",
                                                          workspaceHandler,
                                                          trajectoryHandler,
                                                          do_merged_files_have_correlative_models=True,
                                                          write_frame_number_instead_of_correlative_model_number = keep_frame_number,
                                                          keep_remarks = keep_remarks )

    generatedFiles.append({"description":"Cluster representatives",
                                "path":os.path.abspath(representatives_path),
                                "type":"pdb"})





def save_representatives(representatives,
                         pdb_name,
                         workspace_handler,
                         trajectory_holder,
                         do_merged_files_have_correlative_models,
                         write_frame_number_instead_of_correlative_model_number,
                         keep_remarks = False):
    """
    Saves a pdb file containing the most representative elements of the clustering.

    @param representatives: A list of the representative elements of the clustering we want to extract.

    @param workspace_handler: The workspace handler of this run.

    @param trajectory_holder: The trajectory handler for this run or an array with pdb file paths.

    @param do_merged_files_have_correlative_models: When merging, output file will have models from 0 to M, where M is the total number
    of frames of the merged file.

    @param write_frame_number_instead_of_model_number: When extracting frames, extract those models which number coincides with the
    frame numbers in 'representatives'. Otherwise, extract those models which position coincide with the frame number in
    'representatives'.
    """
    results_directory = workspace_handler["results"]

    # Merge pdbs (in order)
    temporary_merged_trajectory_path = os.path.join(workspace_handler["tmp"],"tmp_merged_trajectory.pdb")

#===========================================================
    # THIS DOES NOT WORK IF USING DCD FILES
#     merge_pdbs(trajectory_holder,
#                temporary_merged_trajectory_path,
#                do_merged_files_have_correlative_models)

    # TEMPORARY HACK TO OVERCOME DCD MERGING BUG

    merged_pdb = trajectory_holder.getMergedStructure()
    prody.writePDB(temporary_merged_trajectory_path, merged_pdb)
#==========================================================

    # Extract frames from the merged pdb
    file_handler_in = open(temporary_merged_trajectory_path,"r")
    file_handler_out = open(os.path.join(results_directory,"%s.pdb"%pdb_name),"w")

    pdb_tools.extract_frames_from_trajectory_sequentially (file_handler_in = file_handler_in,
                                             number_of_frames = pdb_tools.get_number_of_frames(temporary_merged_trajectory_path),
                                             file_handler_out = file_handler_out,
                                             frames_to_save = representatives,
                                             write_frame_number_instead_of_correlative_model_number = write_frame_number_instead_of_correlative_model_number,
                                             keep_header = keep_remarks)
    file_handler_in.close()
    file_handler_out.close()

    return os.path.join(results_directory,"%s.pdb"%pdb_name)

