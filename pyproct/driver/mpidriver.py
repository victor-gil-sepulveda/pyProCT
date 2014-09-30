"""
Created on Mar 26, 2013

@author: victor
"""
from mpi4py import MPI
from pyproct.driver.driver import Driver
from pyproct.driver.workspace.MPIWorkspaceHandler import MPIWorkspaceHandler
from pyproct.driver.time.timerHandler import timed_method
from pyproct.data.MPIDataDriver import MPIDataDriver

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
                self.data_handler, self.matrix_handler = MPIDataDriver.run( self.rank,
                                                                            self.comm,
                                                                            parameters["data"],
                                                                            self.workspaceHandler,
                                                                            Driver.timer,
                                                                            self.generatedFiles)
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
