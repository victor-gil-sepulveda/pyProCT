'''
Created on 07/02/2013

@author: victor
'''
from pyproclust.driver.scheduling.processParallelScheduler import ProcessParallelScheduler
from pyproclust.driver.scheduling.serialScheduler import SerialScheduler


def build_scheduler(scheduler_type, sleep_time, observer, max_processes):
    """
    Factory function for scheduler building.
    @param scheduler_type: Type of the scheduler to be built:
        - "Process/Parallel": A scheduler based on 'multiprocessing'. Can execute tasks in parallel.
        - "Serial": A serial scheduler. Executes its tasks sequentially.
    
    @param sleep_time: The time to wait until a new task is put into the execution queue.
    
    @param observer: The associated observer (can be None)
    
    @param max_processes: Maximum number of simultaneous tasks running.
    
    @return: The scheduler instance.
    """
    if scheduler_type == "Process/Parallel":
        return ProcessParallelScheduler(max_processes, sleep_time, observer)
    elif scheduler_type == "Serial":
        return SerialScheduler(sleep_time, observer)
    else:
        print "[ERROR][ClusteringExplorator::__init__] Not supported scheduler_type ( %s )"%(scheduler_type)
        exit()