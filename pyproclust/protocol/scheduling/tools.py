'''
Created on 07/02/2013

@author: victor
'''
from pyproclust.protocol.scheduling.processParallelScheduler import ProcessParallelScheduler
from pyproclust.protocol.scheduling.serialScheduler import SerialScheduler


def build_scheduler(scheduler_type, sleep_time, observer, max_processes):
    """
    """
    if scheduler_type == "Process/Parallel":
        return ProcessParallelScheduler(max_processes, sleep_time, observer)
    elif scheduler_type == "Serial":
        return SerialScheduler(sleep_time, observer)
    else:
        print "[ERROR][ClusteringExplorator::__init__] Not supported scheduler_type ( %s )"%(scheduler_type)
        exit()