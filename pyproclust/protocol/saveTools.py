'''
Created on 20/09/2012

@author: victor
'''
import shutil
import pyproclust.tools.pdbTools as pdb_tools
import os.path
from pyproclust.driver.handlers.trajectoryHandler import TrajectoryHandler
from pyproclust.tools.pdbTools import get_model_boundaries,\
    repair_MODEL_ENDMDL_tags

def merge_pdbs(trajectories, temporary_merged_trajectory_path, do_merged_files_have_correlative_models):
    pdbs = []
    if trajectories.__class__ == TrajectoryHandler:
        pdbs = [pdb["source"] for pdb in trajectories.pdbs] 
    else:
        pdbs.extend(trajectories)
    
    # Copy the first one (there's at least one)
    shutil.copyfile(pdbs[0], temporary_merged_trajectory_path)
    file_handler_out = open(temporary_merged_trajectory_path, "a")
    for pdb_file in pdbs[1:]:
        # Concat the other files
        file_handler_out.write(open(pdb_file,"r").read())
    file_handler_out.close()
    
    if do_merged_files_have_correlative_models:
        shutil.copyfile(temporary_merged_trajectory_path, 
                        temporary_merged_trajectory_path+".mdl")
        file_handler_in = open(temporary_merged_trajectory_path+".mdl", "r")
        file_handler_out = open(temporary_merged_trajectory_path,"w")
        current_model =  0
        for line in file_handler_in:
            if line[0:5] == "MODEL":
                file_handler_out.write("MODEL "+str(current_model).rjust(8)+"\n")
                current_model = current_model + 1
            else:
                file_handler_out.write(line)
        
        file_handler_out.close()
        file_handler_in.close()
    
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
    merge_pdbs(trajectory_holder, 
               temporary_merged_trajectory_path, 
               do_merged_files_have_correlative_models)
    
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
