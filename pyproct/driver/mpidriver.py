"""
Created on Mar 26, 2013

@author: victor
"""
from mpi4py import MPI
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.driver.driver import Driver
from pyproct.driver.workspace.MPIWorkspaceHandler import MPIWorkspaceHandler
from pyproct.data.dataDriver import DataDriver
from pyproct.driver.time.timerHandler import timed_method
from pyproct.data.matrix.matrixHandler import MatrixHandler

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
                if self.rank == 0:
                    self.data_handler, self.matrix_handler = DataDriver.run(parameters["data"],
                                                                    self.workspaceHandler,
                                                                    Driver.timer,
                                                                    self.generatedFiles)
                else:
                    self.matrix_handler = MatrixHandler(parameters["data"]["matrix"], None)
                self.comm.Barrier()

                matrix_contents = list(self.comm.bcast(self.matrix_handler.distance_matrix.get_data(), root=0))
                if self.rank != 0:
                    self.matrix_handler.distance_matrix = CondensedMatrix(matrix_contents)
                self.comm.Barrier()
                
                if "clustering" in parameters:
                    clustering_results = self.clustering_section(parameters)
                    self.comm.Barrier()

                    if self.rank == 0:
                        self.postprocess(parameters, clustering_results)
                        self.save_results(clustering_results)
                        self.show_summary(parameters, clustering_results)
                        return self.get_best_clustering(clustering_results)
                else:
                    print "[Warning MPIDriver::run] 'clustering' object was not defined in the control script. pyProCT will now stop."
                    self.notify("Driver Finished", "\n"+str(Driver.timer))
            else:
                print "[Warning MPIDriver::run] 'data' object was not defined in the control script. pyProCT will now stop."
                self.notify("MPIDriver Finished", "\n"+str(Driver.timer))

        if self.rank == 0:
            self.notify("MPI-Driver Finished", "\n"+str(self.timer))
