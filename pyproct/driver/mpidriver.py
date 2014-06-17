'''
Created on Mar 26, 2013

@author: victor
'''
from pyproct.driver.handlers.timerHandler import TimerHandler
from pyproct.driver.handlers.trajectoryHandler import TrajectoryHandler
from pyproct.driver.driver import Driver
from mpi4py import MPI
import copy_reg
from pyRMSD.condensedMatrix import CondensedMatrix
import pyRMSD.condensedMatrix
from pyproct.driver.handlers.matrix.matrixHandler import MatrixHandler

class MPIDriver(Driver):

    def __init__(self, observer):
        super(MPIDriver, self).__init__(observer)
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.workspaceHandler = None
        self.trajectoryHandler = None
        self.matrixHandler = None

    def perform_actions(self, parameters):
        best_clustering, clustering_results = self.get_best_clustering(parameters)
        self.comm.Barrier()
        if self.rank == 0:
            self.postprocess(parameters, best_clustering)
            #################################
            # Results are saved to a file
            #################################
            self.save_clustering_results(clustering_results)

            #################################
            # Small summary of the best_cluster is shown
            #################################
            if parameters["clustering"]["generation"]["method"] != "load":
                print "======================="
                print "This is the chosen cluster:"
                print "\t- Used algorithm: ", best_clustering['type']
                print "\t- Number of clusters: ", best_clustering['evaluation']['Number of clusters']
                print "\t- Noise: %.3f %%"%best_clustering['evaluation']['Noise level']

    def run(self, parameters):
        #####################
        # Start timing
        #####################
        self.timer = TimerHandler()
        self.timer.start("Global")

        #####################
        # Workspace Creation
        #####################
        if self.rank == 0:
            self.create_workspace(parameters)
        self.comm.Barrier()
        self.workspaceHandler = self.comm.bcast(self.workspaceHandler, root=0)

        #####################
        # Saving Parameters
        #####################
        if self.rank == 0:
            self.save_parameters_file(parameters)

        #####################
        # Trajectory Loading
        #####################
        self.timer.start("Trajectory Loading")
        self.trajectoryHandler = TrajectoryHandler(parameters, self.observer)
        self.comm.Barrier()
        self.timer.stop("Trajectory Loading")

        ##############################
        # Distance Matrix Generation
        ##############################
        if self.rank == 0:
            self.create_matrix(parameters)
        else:
            self.matrixHandler = MatrixHandler(parameters["data"]["matrix"])
            self.matrixHandler.distance_matrix = CondensedMatrix([0.])
        self.comm.Barrier()
        matrix_contents = list(self.comm.bcast(self.matrixHandler.distance_matrix.get_data(), root=0))
        if self.rank != 0:
            self.matrixHandler.distance_matrix = CondensedMatrix(matrix_contents)

        ##############################
        # Actions
        ##############################
        self.comm.Barrier()
        self.perform_actions(parameters)

        self.timer.stop("Global")
        self.notify("MPI-Driver Finished", "\n"+str(self.timer))

