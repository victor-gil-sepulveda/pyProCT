'''
Created on 20/09/2012

@author: victor
'''
import pickle
import json
import shutil
import pyproclust.tools.pdbTools as pdb_tools
import os.path

def save_representatives(representatives, pdb_name, workspace_handler, trajectory_handler):
    """
    Saves a pdb file containing the most representative elements of the clustering.
    
    @param representatives: A list of the representative frames we want to extract.
    
    @param workspace_handler: The workspace handler of this run.
    
    @param trajectory_handler: The trajectory handler for this run.
    """
    results_directory = workspace_handler["results"]
    temporary_merged_trajectory_path = os.path.join(workspace_handler["tmp"],"tmp_merged_trajectory.pdb")
    pdbs = trajectory_handler.pdbs
    
    # Copy the first one (there's at least one)
    shutil.copyfile(pdbs[0]["source"], temporary_merged_trajectory_path)
    file_handler_in = open(temporary_merged_trajectory_path,"a")
    for pdb_file in pdbs[1:]:
        # Concat the other file
        file_handler_in.write(open(pdb_file["source"],"r").read())
    file_handler_in.close()
    
    # Add 
    file_handler_in = open(temporary_merged_trajectory_path,"r")
    file_handler_out = open(os.path.join(results_directory,"%s.pdb"%pdb_name),"w")

    pdb_tools.extract_frames_from_trajectory(file_handler_in = file_handler_in, 
                                             number_of_frames = pdb_tools.get_number_of_frames(temporary_merged_trajectory_path),
                                             file_handler_out = file_handler_out, 
                                             frames_to_save = representatives,
                                             use_frame_number_as_model = True)
    file_handler_in.close()
    file_handler_out.close()
    
    return os.path.join(results_directory,"%s.pdb"%pdb_name)
  

