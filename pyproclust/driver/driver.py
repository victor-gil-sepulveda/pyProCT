'''
Created on Mar 26, 2013

@author: victor
'''
import os
from pyproclust.driver.handlers.timerHandler import TimerHandler
from pyproclust.driver.handlers.workspaceHandler import WorkspaceHandler
from pyproclust.driver.handlers.trajectoryHandler import TrajectoryHandler
from pyproclust.driver.handlers.matrix.matrixHandler import MatrixHandler
from pyproclust.driver.observer.observable import Observable
import pyproclust.tools.plotTools as plotTools
from pyproclust.protocol.protocol import ClusteringProtocol
import pyproclust.protocol.saveTools as saveTools
from pyproclust.clustering.comparison.distrprob.kullbackLieblerDivergence import KullbackLeiblerDivergence
from pyproclust.driver.compressor.compressor import Compressor
from pyproclust.driver.results.clusteringResultsGatherer import ClusteringResultsGatherer
import json
from pyproclust.clustering.clustering import Clustering

class Driver(Observable):
    
    def __init__(self, observer):
        super(Driver, self).__init__(observer)
    
    def run(self, parameters):
        
        #####################
        # Start timing 
        #####################
        self.timer = TimerHandler()
        self.timer.start("Global")
        
        best_clustering = None
        if parameters["clustering"]["generation"]["method"] == "generate":
            #####################
            # Workspace Creation 
            #####################
            self.workspaceHandler = WorkspaceHandler(parameters["workspace"], self.observer)
            self.workspaceHandler.create_directories()
            
            #####################
            # Saving Parameters 
            #####################
            parameters_file_path = os.path.join(self.workspaceHandler["results"],"parameters.json")
            open(parameters_file_path,"w").write(json.dumps(parameters.params_dic,
                                                              sort_keys=False,
                                                              indent=4,
                                                              separators=(',', ': ')))
            
            self.generatedFiles = [{"description":"Parameters file", "path":parameters_file_path,"type":"text"}]
            
            #####################
            # Trajectory Loading
            #####################
            self.timer.start("Trajectory Loading")
            self.trajectoryHandler = TrajectoryHandler(parameters["global"], parameters["matrix"], self.observer)
            self.timer.stop("Trajectory Loading")
            
            ##############################
            # Distance Matrix Generation
            ##############################
            self.matrixHandler = MatrixHandler(parameters["matrix"])
            self.notify("Matrix calculation",[])
            self.timer.start("Matrix Generation")
            self.matrixHandler.create_matrix(self.trajectoryHandler)
            statistics_file_path = self.matrixHandler.save_statistics(self.workspaceHandler["matrix"])
            self.generatedFiles.append({"description":"Matrix statistics", "path":statistics_file_path,"type":"text"})
            self.timer.stop("Matrix Generation")
            self.timer.start("Matrix Save")
            self.matrixHandler.save_matrix(os.path.join(self.workspaceHandler["matrix"],parameters["matrix"]["filename"]))
            self.timer.stop("Matrix Save")
    
            #########################
            # Matrix plot
            #########################
            if "image" in parameters["matrix"].keys() :
                self.timer.start("Matrix Imaging")
                matrix_image_file_path = os.path.join(self.workspaceHandler["matrix"], parameters["matrix"]["image"]["filename"])
                plotTools.matrixToImage(  self.matrixHandler.distance_matrix,
                                          matrix_image_file_path,
                                          max_dim = parameters["matrix"]["image"]["dimension"],
                                          observer = self.observer)
                self.generatedFiles.append({"description":"Matrix image", "path":matrix_image_file_path,"type":"image"})
                self.timer.stop("Matrix Imaging")
             
            ##############################
            # Do the actual clustering
            ##############################
            print "performing protocol"
            clustering_results = ClusteringProtocol(self.timer, self.observer).run(parameters,
                                                                                self.matrixHandler,
                                                                                self.workspaceHandler,
                                                                                self.trajectoryHandler)
            best_clustering_id, selected, not_selected, scores = (None, {},{},{})
            best_clustering = None
            
            if clustering_results is not None:
                best_clustering_id, selected, not_selected, scores = clustering_results
                best_clustering = selected[best_clustering_id]
            else:
                pass
                #SHUTDOWN, NO CLUSTER
        
        if parameters["clustering"]["generation"]["method"] == "load":
            best_clustering = Clustering.from_dic(parameters["clustering"]["generation"]);
              
        ##############################
        # Specialized post-processing
        ##############################
        action_type = parameters["global"]["action"]["type"]
        if action_type == "clustering" or action_type == "advanced":
            self.timer.stop("Global")

            #################################
            # Abort if no clusters were found
            #################################
            if best_clustering is None:
                self.notify("SHUTDOWN", "The clustering search found no clusterings. Relax evaluation constraints.")
                print "[FATAL ClusteringProtocol::run] The clustering search found no clusterings. Exiting..."
                exit()
                
            ##############################
            # Saving representatives
            ##############################
            medoids = best_clustering["clustering"].get_medoids(self.matrixHandler.distance_matrix)
            
            # Set prototypes and ids (medoids are ordered)
            for i in range(len(best_clustering["clustering"].clusters)):
                best_clustering["clustering"].clusters[i].prototype = medoids[i]
            
            representatives_path = saveTools.save_representatives(medoids, 
                                                                   "representatives",
                                                                   self.workspaceHandler, 
                                                                   self.trajectoryHandler,
                                                                   do_merged_files_have_correlative_models = True,
                                                                   write_frame_number_instead_of_correlative_model_number = True)
            
            self.generatedFiles.append({"description":"Cluster central conformations", "path":representatives_path, "type":"pdb"})

            #################################
            # Results are saved to a file
            #################################
            results_path = os.path.join(self.workspaceHandler["results"],"results.json")
            self.generatedFiles.append({"description":"Results file", "path":results_path, "type":"text"})

            json_results = ClusteringResultsGatherer().gather(self.timer, 
                                                              self.trajectoryHandler, 
                                                              self.workspaceHandler,
                                                              clustering_results,
                                                              self.generatedFiles)
            # Results are first added and saved later to avoid metareferences :D
            open(results_path,"w").write(json_results)
            
        elif action_type == "comparison":
            ############################################
            # Distribution analysis
            ############################################
            self.timer.start("KL divergence")
            klDiv = KullbackLeiblerDivergence(self.trajectoryHandler.pdbs, self.matrixHandler.distance_matrix)
            kl_file_path = os.path.join(self.workspaceHandler["matrix"],"kullback_liebler_divergence")
            klDiv.save(kl_file_path)
            self.generatedFiles.append({"description":"Kullback-Leibler divergence", "path":matrix_image_file_path,"type":"text"})
            
            self.timer.stop("KL divergence")
            self.timer.stop("Global")
        
        elif action_type == "compression":
            ############################################
            # Compress
            ############################################
            self.timer.start("Compression")
            compressor = Compressor(parameters["global"]["action"]["parameters"])
            compressed_file_path = compressor.compress(best_clustering["clustering"], "compressed_pdb",
                               self.workspaceHandler, 
                               self.trajectoryHandler,
                               self.matrixHandler)
            self.generatedFiles.append({"description":"Compressed file", "path":compressed_file_path,"type":"pdb"})
            self.timer.stop("Compression")
            self.timer.stop("Global")
        
        print self.timer
        
