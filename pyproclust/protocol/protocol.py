'''
Created on 21/05/2012

@author: victor
'''
import pyproclust.driver.scheduling.tools as scheduling_tools
from pyproclust.protocol.exploration.clusteringExplorator import ClusteringExplorator
from pyproclust.clustering.filtering.clusteringFilter import ClusteringFilter
from pyproclust.clustering.analysis.picklingAnalysisRunner import PicklingAnalysisRunner
from pyproclust.clustering.analysis.analysisPopulator import AnalysisPopulator
from pyproclust.clustering.selection.bestClusteringSelector import BestClusteringSelector
import pyproclust.protocol.saveTools as save_tools
from pyproclust.driver.observer.observable import Observable

class ClusteringProtocol(Observable):

    def __init__(self, observer):
        super(ClusteringProtocol, self).__init__(observer)
        
    def run(self, parameters):
        ############################
        # Clustering exploration
        ############################
        self.timer.start("Clustering Exploration")
        clusterings  = ClusteringExplorator(parameters, 
                                           self.matrixHandler, 
                                           self.workspaceHandler, 
                                           scheduling_tools.build_scheduler("Process/Parallel",
                                                           parameters["control"]["algorithm_scheduler_sleep_time"],
                                                           self.observer,
                                                           parameters["control"]["number_of_processors"]), 
                                           self.observer).run()
        self.notify("Clusterings Created", {"number_of_clusters":len(clusterings)})
        self.timer.stop("Clustering Exploration")
        
        ######################
        # First filtering
        ######################
        self.timer.start("Clustering Filtering")
        selected_clusterings, not_selected_clusterings = ClusteringFilter(parameters["evaluation"], 
                                                                          self.matrixHandler).filter(clusterings)
            
        self.notify("Filter", {"selected":selected_clusterings,"not_selected":not_selected_clusterings})
        self.timer.stop("Clustering Filtering")
        
        if selected_clusterings == []:
            self.notify("SHUTDOWN", "The clustering search gave no clusterings. Relax evaluation constraints.")
        else:      
            ######################
            # Clustering scoring
            ######################
            self.timer.start("Evaluation")
            analyzer = PicklingAnalysisRunner(scheduling_tools.build_scheduler("Process/Parallel",
                                                           parameters["control"]["algorithm_scheduler_sleep_time"],
                                                           self.observer,
                                                           parameters["control"]["number_of_processors"]),
                                              parameters,
                                              selected_clusterings,
                                              AnalysisPopulator(self.matrixHandler, 
                                                                self.trajectoryHandler,
                                                                parameters))
            
            analyzer.evaluate()
            self.timer.stop("Evaluation")
            
            ######################
            # Choose the best clustering
            ######################
            self.timer.start("Selection")
            best_clustering_id, best_criteria_id, all_scores = BestClusteringSelector(parameters).choose_best(selected_clusterings)
            self.notify("BEST", (best_clustering_id, best_criteria_id, all_scores[best_criteria_id][best_clustering_id]))
            self.timer.stop("Selection")

            ######################
            # Save results
            ######################
            save_tools.save_cluster_info(self.workspaceHandler["results"]+"/selected", selected_clusterings)
            save_tools.save_cluster_info(self.workspaceHandler["results"]+"/not_selected", not_selected_clusterings)
            save_tools.save_best_clusters_and_scores(best_clustering_id, 
                                                     best_criteria_id, 
                                                     all_scores, 
                                                     self.workspaceHandler["results"]+"/scores")
            
            
            
#             #########################
#             # Get statistics etc...
#             #########################
#             if protocol_params.most_representative_pdb_file != "":
#                 save_most_representative(protocol_params,best_clustering,\
#                                         self.matrixHandler.distance_matrix,\
#                                         self.workspaceHandler.tmp_path,\
#                                         self.workspaceHandler.results_path)
#             
#             
#             analyzer = ClusteringStatisticalAnalyzer(best_clustering,\
#                                                      protocol_params.pdb1,\
#                                                      protocol_params.pdb2,\
#                                                      self.matrixHandler.distance_matrix,\
#                                                      protocol_params.shallWeCompareTrajectories())
#             analyzer.per_cluster_analytics()
#             
#             analyzer.per_clustering_analytics()
#             
#             plotGenerator = ClusteringPlotsGenerator(analyzer)
#             
#             big_plot = plotGenerator.generate_and_compose_big_plot(composition_size= (1024,720), max_radius = 150, ball_horizontal_separation = 50, ball_vertical_separation = 100)
#             
#             if protocol_params.shallWeCompareTrajectories():
#                 small_plot = plotGenerator.generate_and_compose_small_plot(composition_size= (400,700))
#             else:
#                 small_plot = plotGenerator.generate_and_compose_small_plot(composition_size= (400,300))
#             
#             big_plot.save(self.workspaceHandler.results_path+"/analysis_plot_big.png")
#             small_plot.save(self.workspaceHandler.results_path+"/analysis_plot_small.png")
#             
#         global_time_end = time.time()
#         self.htmlReport.report["Timing"] += 'All the process took %0.3f s\n' % (global_time_end-global_time_start)
#         
#         time_file = open(self.workspaceHandler.tmp_path+"/timing.txt","w")
#         time_file.write(self.htmlReport.report["Timing"])
#         time_file.close()
#         
#         
#         #HTML REPORT
#         html_file = open(self.workspaceHandler.results_path+"/report.html","w")
#         html_file.write(self.htmlReport.generateHTML())
#         html_file.close()
#         
#         #IMAGES FOR REPORT
#         self.htmlReport.report["Image Paths"]["kl"] = self.workspaceHandler.matrix_path+"/rmsd_distrib.png"
#         self.htmlReport.report["Image Paths"]["matrix"] = self.workspaceHandler.matrix_path+"/matrix_plot.png"
#         self.htmlReport.report["Image Paths"]["clustering_small"] = self.workspaceHandler.results_path+"/analysis_plot_small.png"
#         self.htmlReport.report["Image Paths"]["clustering_big"] = self.workspaceHandler.results_path+"/analysis_plot_big.png"
#         self.htmlReport.create_thumbnails(self.workspaceHandler.results_path)
#     
#    
        
    