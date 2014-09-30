"""
Created on 30/9/2014

@author: victor
"""
from pyproct.data.dataDriver import DataDriver
from pyproct.data.matrix.matrixHandler import MatrixHandler
from pyRMSD.condensedMatrix import CondensedMatrix
import os

class MPIDataDriver(DataDriver): 
    def __init__(self):
        super(MPIDataDriver, self).__init__()
    
    @classmethod
    def run(cls, rank, mpi_comm, parameters, workspace_handler, timer, generated_files):
        DataDriver.timer = timer
        
        # Wait until the trajectories are loaded. 
        # Ideally this and a lot of other things may be done
        # in parallel too.
        # TODO: Atomic looks not to be pickable, which doesn't allow
        # us to directly share the data_handler
        data_handler = cls.load_data(parameters)
        mpi_comm.Barrier() 

        # Load the trajectories        
        if rank == 0:
            matrix_handler = cls.calc_matrix(data_handler, 
                                         parameters["matrix"])
            matrix_contents = matrix_handler.distance_matrix.get_data()
        else:
            matrix_contents = None
            matrix_handler = MatrixHandler(None, parameters["matrix"])
        mpi_comm.Barrier() 
        # Wait until matrix is calculated.
        
        # Then broadcast its contents
        matrix_contents = mpi_comm.bcast(matrix_contents, root=0)
        if rank != 0:
            matrix_handler.distance_matrix = CondensedMatrix(matrix_contents)
        
        # We expect to have at least 2 processes running (if not, the process will freeze)
        if rank == 0:
            # Save statistics
            statistics_file_path = matrix_handler.save_statistics(workspace_handler["matrix"])
            generated_files.append({
                                    "description":"Matrix statistics",
                                    "path":os.path.abspath(statistics_file_path),
                                    "type":"text"
            })
    
            # Save matrix contents
            if "filename" in parameters["matrix"]:
                cls.save_matrix(matrix_handler, 
                                workspace_handler,
                                parameters["matrix"])
        if rank == 1:
            # Plot matrix
            if "image" in parameters["matrix"]:
                cls.plot_matrix(matrix_handler, 
                                workspace_handler,
                                parameters["matrix"], 
                                generated_files)
        
        return data_handler, matrix_handler 