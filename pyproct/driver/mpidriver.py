"""
Created on Mar 26, 2013

@author: victor
"""
from mpi4py import MPI
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.tools.commonTools import timed_method
from pyproct.driver.driver import Driver
from pyproct.driver.workspace.MPIWorkspaceHandler import MPIWorkspaceHandler
from pyproct.data.proteins.matrix.matrixHandler import MatrixHandler

class MPIDriver(Driver):
    """
    MPI version of the driver.
    """
    def __init__(self, observer):
        super(MPIDriver, self).__init__(observer)
        self.comm = MPI.COMM_WORLD
        self.nprocs = self.comm.Get_size()
        self.rank = self.comm.Get_rank()
        self.workspaceHandler = None
        self.trajectoryHandler = None
        self.matrixHandler = None

    @timed_method(Driver.timer, "Global")
    def run(self, parameters):

        with MPIWorkspaceHandler(self.rank, 0, parameters["global"]["workspace"], self.observer) as self.workspaceHandler:
            self.comm.Barrier()

            if self.rank == 0:
                self.save_parameters_file(parameters)

            if "data" in parameters:
                self.data_section(parameters)

                if "clustering" in parameters:
                    clustering_results = self.clustering_section(parameters)
                    self.comm.Barrier()

                    if self.rank == 0:
                        self.postprocess(parameters, clustering_results)
                        self.save_results(clustering_results)
                        self.show_summary(parameters, clustering_results)

        if self.rank == 0:
            self.notify("MPI-Driver Finished", "\n"+str(self.timer))

    def data_section(self, parameters):
        """
        MPI implementation of the data protocol.
        """
        matrix_parameters = parameters["data"]["matrix"]

        self.load_trajectory(parameters)
        self.comm.Barrier()

        if self.rank == 0:
            self.calculate_matrix(parameters)
        else:
            self.matrixHandler = MatrixHandler(matrix_parameters)
            self.matrixHandler.distance_matrix = CondensedMatrix([0.]) # fake matrix
        self.comm.Barrier()

        matrix_contents = list(self.comm.bcast(self.matrixHandler.distance_matrix.get_data(), root=0))
        if self.rank != 0:
            self.matrixHandler.distance_matrix = CondensedMatrix(matrix_contents)

        if self.nprocs > 1:

            if self.rank == 0:
                if "filename" in matrix_parameters:
                    self.save_matrix(parameters)

            if self.rank == 1:
                if "image" in matrix_parameters:
                    self.plot_matrix(parameters)
        self.comm.Barrier()



