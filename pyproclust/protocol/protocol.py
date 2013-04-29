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
import os
from pyproclust.protocol.exploration.algorithmRunParametersGenerator import AlgorithmRunParametersGenerator

class ClusteringProtocol(Observable):

    def __init__(self, timer, observer):
        super(ClusteringProtocol, self).__init__(observer)
        self.timer = timer
        
    def run(self, clustering_parameters, matrixHandler, workspaceHandler, trajectoryHandler):
        
        ############################
        # Clustering exploration
        ############################
        self.timer.start("Clustering Exploration")
        clusterings  = ClusteringExplorator(
                                            clustering_parameters,
                                            matrixHandler,
                                            workspaceHandler,
                                            scheduling_tools.build_scheduler(
                                                                           "Process/Parallel",
                                                                           clustering_parameters["clustering"]["control"]["algorithm_scheduler_sleep_time"],
                                                                           self.observer,
                                                                           clustering_parameters["clustering"]["control"]["number_of_processors"]), 
                                            AlgorithmRunParametersGenerator(
                                                                            clustering_parameters,
                                                                            matrixHandler),
                                            self.observer).run()
                                            
        self.notify("Clusterings Created", {"number_of_clusters":len(clusterings)})
        self.timer.stop("Clustering Exploration")
        
        ######################
        # First filtering
        ######################
        self.timer.start("Clustering Filtering")
        selected_clusterings, not_selected_clusterings = ClusteringFilter(clustering_parameters["evaluation"], 
                                                                          matrixHandler).filter(clusterings)
            
        self.notify("Filter", {"selected":selected_clusterings,"not_selected":not_selected_clusterings})
        self.timer.stop("Clustering Filtering")
        
        if selected_clusterings == {}:
            return None
        print "SELECTED", selected_clusterings
        ######################
        # Clustering scoring
        ######################
        self.timer.start("Evaluation")
        analyzer = PicklingAnalysisRunner(scheduling_tools.build_scheduler("Process/Parallel",
                                                       clustering_parameters["clustering"]["control"]["evaluation_scheduler_sleep_time"],
                                                       self.observer,
                                                       clustering_parameters["clustering"]["control"]["number_of_processors"]),
                                          selected_clusterings,
                                          AnalysisPopulator(matrixHandler, 
                                                            trajectoryHandler,
                                                            clustering_parameters))
        
        analyzer.evaluate()
        self.timer.stop("Evaluation")
        
        ######################
        # Choose the best clustering
        ######################
        self.timer.start("Selection")
        best_clustering_id, best_criteria_id, all_scores = BestClusteringSelector(clustering_parameters).choose_best(selected_clusterings)
        self.timer.stop("Selection")

        print "BC",best_clustering_id, best_criteria_id, all_scores
        ######################
        # Save results
        ######################
        save_tools.save_cluster_info(os.path.join(workspaceHandler["results"],"selected"), selected_clusterings)
        save_tools.save_cluster_info(os.path.join(workspaceHandler["results"],"not_selected") , not_selected_clusterings)
        
        save_tools.save_best_clusters_and_scores(best_clustering_id, 
                                                 best_criteria_id, 
                                                 all_scores, 
                                                 os.path.join(workspaceHandler["results"],"scores"))
        
        return selected_clusterings[best_clustering_id]
            
#             
#             analyzer = ClusteringStatisticalAnalyzer(best_clustering,\
#                                                      protocol_params.pdb1,\
#                                                      protocol_params.pdb2,\
#                                                      matrixHandler.distance_matrix,\
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
#             big_plot.save(workspaceHandler.results_path+"/analysis_plot_big.png")
#             small_plot.save(workspaceHandler.results_path+"/analysis_plot_small.png")
#             
#         global_time_end = time.time()
#         self.htmlReport.report["Timing"] += 'All the process took %0.3f s\n' % (global_time_end-global_time_start)
#         
#         time_file = open(workspaceHandler.tmp_path+"/timing.txt","w")
#         time_file.write(self.htmlReport.report["Timing"])
#         time_file.close()
#         
#         
#         #HTML REPORT
#         html_file = open(workspaceHandler.results_path+"/report.html","w")
#         html_file.write(self.htmlReport.generateHTML())
#         html_file.close()
#         
#         #IMAGES FOR REPORT
#         self.htmlReport.report["Image Paths"]["kl"] = workspaceHandler.matrix_path+"/rmsd_distrib.png"
#         self.htmlReport.report["Image Paths"]["matrix"] = workspaceHandler.matrix_path+"/matrix_plot.png"
#         self.htmlReport.report["Image Paths"]["clustering_small"] = workspaceHandler.results_path+"/analysis_plot_small.png"
#         self.htmlReport.report["Image Paths"]["clustering_big"] = workspaceHandler.results_path+"/analysis_plot_big.png"
#         self.htmlReport.create_thumbnails(workspaceHandler.results_path)
#     
#    
        
    