'''
Created on Mar 26, 2013

@author: victor
'''
from pyproclust.protocol.handlers.timerHandler import TimerHandler
from pyproclust.protocol.handlers.workspaceHandler import WorkspaceHandler
from pyproclust.protocol.handlers.trajectoryHandler import TrajectoryHandler
from pyproclust.protocol.handlers.matrix.matrixHandler import MatrixHandler
from pyproclust.clustering.comparison.distrprob.kullbackLieblerDivergence import KullbackLeiblerDivergence
import pyproclust.tools.plotTools as plotTools
from pyproclust.protocol.observer.observable import Observable

class driver(Observable):
    '''
    classdocs
    '''


    def __init__(self, observer):
        super(driver, self).__init__(observer)
    
    def run(self,parameters):
        #####################
        # Start timing 
        #####################
        self.timer = TimerHandler()
        self.timer.start("Global")
        
        #####################
        # Create workspace 
        #####################
        self.workspaceHandler = WorkspaceHandler(parameters, self.observer)
        self.workspaceHandler.create_directories()

        #####################
        # Loading trajectory 
        #####################
        self.timer.start("Trajectory Loading")
        self.trajectoryHandler = TrajectoryHandler(parameters, self.observer)
        self.timer.stop("Trajectory Loading")
        
        ##############################
        # Obtaining the distance matrix
        ##############################
        self.matrixHandler = MatrixHandler(parameters["matrix"], parameters)
        
        self.timer.start("Matrix Generation")
        self.matrixHandler.create_matrix()
        self.matrixHandler.save_statistics(parameters["matrix"]["matrix_path"])
        self.timer.stop("Matrix Generation")
           
        if parameters["matrix"]["save_matrix"]:
            self.timer.start("Matrix Save")
            self.matrixHandler.save_matrix(parameters["matrix"]["store_matrix_path"])
            self.timer.stop("Matrix Save")
        
#        ############################################
#        # Distribution analysis
#        ############################################
#        if parameters["matrix"]["action"] == "comparison":
#            self.timer.start("KL divergence")
#            klDiv = KullbackLeiblerDivergence(self.trajectoryHandler.pdbs, self.matrixHandler.distance_matrix)
#            klDiv.save(self.workspaceHandler["matrix"]+"/kullback_liebler_divergence")
#            self.timer.stop("KL divergence")
#        print self.timer.get_elapsed()
#        
#        #########################
#        # Matrix plot
#        #########################
#        self.timer.start("Matrix Imaging")
#        plotTools.matrixToImage(  self.matrixHandler.distance_matrix,
#                                  self.workspaceHandler["matrix"]+"/matrix_plot.png",
#                                  max_dim = 1000,
#                                  observer = self.observer)
#        self.timer.stop("Matrix Imaging")
        
        
        