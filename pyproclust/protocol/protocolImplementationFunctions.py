'''
Created on 20/09/2012

@author: victor
'''
import pyproclust.tools.commonTools as common
import pyproclust.tools.pdbTools as pdb_tools
import pickle
from pyproclust.clustering.filtering.clusteringFilter import ClusteringFilter
from pyproclust.clustering.analysis.picklingParallelAnalysisRunner import PicklingParallelAnalysisRunner
from pyproclust.clustering.analysis.analysisPopulator import AnalysisPopulator

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

  
def save_results(protocol_params,directory,string_results,results_pack):
    if (protocol_params.report_file != ""):
        reports_file_handler = open(directory+"/"+protocol_params.report_file+".txt","w")
        reports_file_handler.write(string_results)
        reports_file_handler.close()
        result_pack_file_handler = open(directory+"/"+protocol_params.report_file+".bin","w")
        pickle.dump(results_pack,result_pack_file_handler)
        result_pack_file_handler.close()
    else:
        print string_results
        
def clustering_scoring(clustering_info_structures, parameters, condensed_matrix, pdb_structure):
    
    analyzer = PicklingParallelAnalysisRunner(scheduling_tools.build_scheduler("Process/Parallel",
                                                           parameters["control"]["algorithm_scheduler_sleep_time"],
                                                           self.observer,
                                                           parameters["control"]["number_of_processors"])
    
    AnalysisPopulator(analyzer, condensed_matrix, pdb_structure,protocol_params.evaluation_types)
    for c in already_filtered_clusterings:
        analyzer.run_analysis_for(c)
    return  analyzer.generate_report()