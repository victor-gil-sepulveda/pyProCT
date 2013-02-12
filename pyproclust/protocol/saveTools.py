'''
Created on 20/09/2012

@author: victor
'''
import pickle
import json
import shutil
import pyproclust.tools.pdbTools as pdb_tools

def save_cluster_info(filename, cluster_info):
    """
    Pickles the cluster_info structure which (depending at which point of the protocol this function is called) may have
    all needed info.
    """
    file_handler = open(filename+".bin", "w")
    pickle.dump(cluster_info, file_handler)
    file_handler.close()
    
def save_best_clusters_and_scores(best_clustering_id, best_criteria_id, all_scores, filename):
    """
    Gathers all the important scoring features and writes them to disk.
    """
    scoring_dic = {
                  "best_cluster": best_clustering_id,
                  "best_criteria": best_criteria_id,
                  "best_score": all_scores[best_criteria_id][best_clustering_id],
                  "scores": all_scores
    }
    open(filename+".json","w").write(json.dumps(scoring_dic, sort_keys=True, indent=4, separators=(',', ': ')))

def save_representatives(representatives, distance_matrix, workspace_handler, trajectory_handler):
    """
    Saves a pdb file containing the most representative elements of the clustering.
    
    @param representatives: A list of the representative frames we want to extract.
    
    @param distance_matrix: The distance matrix used to get the clustering.
    
    @param workspace_handler: The workspace handler of this run.
    
    @param trajectory_handler: The trajectory handler for this run.
    """
    results_directory = workspace_handler["results"]
    temporary_merged_trajectory_path = workspace_handler["tmp"]+"/tmp_merged_trajectory.pdb"
    pdbs = trajectory_handler.pdbs
    
    # Copy the first one (there's at least one)
    shutil.copyfile(pdbs[0], temporary_merged_trajectory_path)
    file_handler_in = open(temporary_merged_trajectory_path,"a")
    for pdb_file in pdbs[1:]:
        # Concat the other file
        file_handler_in.write(open(pdb_file,"r").read())
    file_handler_in.close()
    
    # Add 
    file_handler_in = open(temporary_merged_trajectory_path,"r")
    file_handler_out = open(results_directory+"/representatives.pdb","w")
    pdb_tools.extract_frames_from_trajectory(file_handler_in, distance_matrix.row_length, file_handler_out, representatives)
    file_handler_in.close()
    file_handler_out.close()

  

