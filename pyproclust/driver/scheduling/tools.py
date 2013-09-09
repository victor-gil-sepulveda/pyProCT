'''
Created on 07/02/2013

@author: victor
'''
from pyscheduler.serialScheduler import SerialScheduler
from pyscheduler.processParallelScheduler import ProcessParallelScheduler

def build_scheduler(scheduler_type, observer, max_processes):
    """
    Factory function for scheduler building.
    @param scheduler_type: Type of the scheduler to be built:
        - "Process/Parallel": A scheduler based on 'multiprocessing'. Can execute tasks in parallel.
        - "MPI/Parallel": This one uses MPI to execute tasks in parallel. Only usable with an MPI driver.
        - "Serial": A serial scheduler. Executes its tasks sequentially.
    
    @param sleep_time: The time to wait until a new task is put into the execution queue.
    
    @param observer: The associated observer (can be None)
    
    @param max_processes: Maximum number of simultaneous tasks running.
    
    @return: The scheduler instance.
    """
    if scheduler_type == "Process/Parallel":
        return ProcessParallelScheduler(max_processes)
    elif scheduler_type == "MPI/Parallel":
        from pyscheduler.MPIParallelScheduler import MPIParallelScheduler # to avoid unneeded call to mpi_init
        return MPIParallelScheduler(share_results_with_all_processes=True)
    elif scheduler_type == "Serial":
        return SerialScheduler()
    else:
        print "[ERROR][ClusteringExplorator::__init__] Not supported scheduler_type ( %s )"%(scheduler_type)
        exit()