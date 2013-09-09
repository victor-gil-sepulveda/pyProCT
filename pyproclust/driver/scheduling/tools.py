'''
Created on 07/02/2013

@author: victor
'''
from pyscheduler.serialScheduler import SerialScheduler
from pyscheduler.processParallelScheduler import ProcessParallelScheduler
import time, datetime

def send_message_to_observer( observer, tag, task_name):
    timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('[%H:%M:%S]')
    if task_name is not None:
        observer.notify("Scheduler", tag, task_name+" "+timestamp)
    else:
        observer.notify("Scheduler", tag, timestamp)

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
    
    # Define functions
    external_functions = {
                          'task_started':{
                                        'function':send_message_to_observer,
                                        'kwargs':{
                                                  'observer':observer,
                                                  'tag':'Task Started'
                                                  }
                                        },
                          'task_ended':{
                                      'function':send_message_to_observer,
                                      'kwargs':{
                                                'observer':observer,
                                                'tag':'Task Ended'
                                                }
                                      },
                          'scheduling_ended':{
                                             'function':send_message_to_observer,
                                             'kwargs':{
                                                       'observer':observer,
                                                       'tag':'Scheduler Ended'
                                                       }
                                             }
                          }
    
    if scheduler_type == "Process/Parallel":
        return ProcessParallelScheduler(max_processes,external_functions)
    elif scheduler_type == "MPI/Parallel":
        from pyscheduler.MPIParallelScheduler import MPIParallelScheduler # to avoid unneeded call to mpi_init
        return MPIParallelScheduler(share_results_with_all_processes=True, functions = external_functions)
    elif scheduler_type == "Serial":
        return SerialScheduler(external_functions)
    else:
        print "[ERROR][ClusteringExplorator::__init__] Not supported scheduler_type ( %s )"%(scheduler_type)
        exit()