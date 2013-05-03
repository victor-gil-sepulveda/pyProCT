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

class Driver(Observable):
    def __init__(self, observer):
        super(Driver, self).__init__(observer)
    
    def run(self,parameters):
        #####################
        # Start timing 
        #####################
        self.timer = TimerHandler()
        self.timer.start("Global")
        
        #####################
        # Workspace Creation 
        #####################
        self.workspaceHandler = WorkspaceHandler(parameters["workspace"], self.observer)
        self.workspaceHandler.create_directories()

        #####################
        # Trajectory Loading
        #####################
        self.timer.start("Trajectory Loading")
        self.trajectoryHandler = TrajectoryHandler(parameters["global"], self.observer)
        self.timer.stop("Trajectory Loading")
        
        ##############################
        # Distance Matrix Generation
        ##############################
        self.matrixHandler = MatrixHandler(parameters["matrix"])
        self.notify("Matrix calculation",[])
        self.timer.start("Matrix Generation")
        self.matrixHandler.create_matrix(self.trajectoryHandler)
        self.matrixHandler.save_statistics( self.workspaceHandler["matrix"])
        self.timer.stop("Matrix Generation")
        self.timer.start("Matrix Save")
        self.matrixHandler.save_matrix(os.path.join(self.workspaceHandler["matrix"],parameters["matrix"]["filename"]))
        self.timer.stop("Matrix Save")


        action_type = parameters["global"]["action"]["type"]
        
        #########################
        # Matrix plot
        #########################
        if "image" in parameters["matrix"].keys() :
            self.timer.start("Matrix Imaging")
            plotTools.matrixToImage(  self.matrixHandler.distance_matrix,
                                      os.path.join(self.workspaceHandler["matrix"], parameters["matrix"]["image"]["filename"]),
                                      max_dim = parameters["matrix"]["image"]["dimension"],
                                      observer = self.observer)
               
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
        
        if clustering_results[0] is not None:
            best_clustering_id, selected, not_selected, scores = clustering_results
            best_clustering = selected[best_clustering_id]
            
        ##############################
        # Specialized post-processing
        ##############################
        if action_type == "clustering":
            #################################
            # Saving results
            #################################
            json_results = ClusteringResultsGatherer().gather(self.timer, self.trajectoryHandler, clustering_results)
            open(os.path.join(self.workspaceHandler["results"],"results.json"),"w").write(json_results)
            
            #################################
            # Abort if no clusters found
            #################################
            if best_clustering is None:
                self.notify("SHUTDOWN", "The clustering search found no clusterings. Relax evaluation constraints.")
                print "[FATAL ClusteringProtocol::run] The clustering search found no clusterings. Exiting..."
                exit()
                
            ##############################
            # Saving representatives
            ##############################
            saveTools.save_representatives(best_clustering["clustering"].get_medoids(self.matrixHandler.distance_matrix), 
                                           "representatives",
                                           self.workspaceHandler, 
                                           self.trajectoryHandler)
            self.timer.stop("Global")
            
        elif action_type == "comparison":
            ############################################
            # Distribution analysis
            ############################################
            self.timer.start("KL divergence")
            klDiv = KullbackLeiblerDivergence(self.trajectoryHandler.pdbs, self.matrixHandler.distance_matrix)
            klDiv.save(os.path.join(self.workspaceHandler["matrix"],"kullback_liebler_divergence"))
            self.timer.stop("KL divergence")
            self.timer.stop("Global")
        
        elif action_type == "compression":
            ############################################
            # Compress
            ############################################
            self.timer.start("Compression")
            compressor = Compressor(parameters["global"]["action"]["parameters"])
            compressor.compress(best_clustering["clustering"], "compressed_pdb",
                               self.workspaceHandler, 
                               self.trajectoryHandler,
                               self.matrixHandler)
            self.timer.stop("Compression")
            self.timer.stop("Global")
        
        print self.timer
        
