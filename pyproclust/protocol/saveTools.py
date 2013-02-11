'''
Created on 20/09/2012

@author: victor
'''
import pyproclust.tools.commonTools as common
import pyproclust.tools.pdbTools as pdb_tools
import pickle
import json


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
    Gathers all the important scoring 
    """
    scoring_dic = {
                  "best_cluster": best_clustering_id,
                  "best_criteria": best_criteria_id,
                  "best_score": all_scores[best_criteria_id][best_clustering_id],
                  "scores": all_scores
    }
    open(filename+".json","w").write(json.dumps(scoring_dic, sort_keys=True, indent=4, separators=(',', ': ')))

def save_most_representative_with_statistical_significance():
    pass

# get_medoids(self, distance_matrix)
# get_proportional_size_representatives(self, number_of_structures, distance_matrix)

def save_most_representative(protocol_params, clustering, distance_matrix, tmp_directory,results_directory):
    medoids= []
    for c in clustering.clusters:
        medoids.append(c.calculate_medoid(distance_matrix))
    
    temporary_merged_trajectory_path = tmp_directory+"/tmp_merged_trajectory.pdb"
    file_handler_in = None
    if protocol_params.shallWeMergetrajectories(): 
        # We want to merge both trajectories and filter them (we're only interested in alphas right now)
        common.print_and_flush("Merging trajectories ...")
        file_handler_in = open(temporary_merged_trajectory_path,"w")
        pdb1_fh = open(protocol_params.pdb1)
        pdb2_fh = open(protocol_params.pdb2)
        common.merge_files([pdb1_fh,pdb2_fh],file_handler_in,False)
        file_handler_in.close()
        pdb1_fh.close()
        pdb2_fh.close()
        common.print_and_flush(" Done\n")
        file_handler_in = open(temporary_merged_trajectory_path,"r")
    else:
        # We are analyzing only one trajectory
        file_handler_in = open(protocol_params.pdb1)
    
    file_handler_out = open(results_directory+"/"+protocol_params.most_representative_pdb_file,"w")
    pdb_tools.extract_frames_from_trajectory(file_handler_in, distance_matrix.row_length, file_handler_out, medoids)
    file_handler_in.close()
    file_handler_out.close()

  

