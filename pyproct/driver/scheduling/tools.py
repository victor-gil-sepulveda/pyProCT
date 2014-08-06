"""
Created on 07/02/2013

@author: victor
"""
from pyscheduler.serialScheduler import SerialScheduler
from pyscheduler.processParallelScheduler import ProcessParallelScheduler
import time, datetime
import multiprocessing

def send_message_to_observer( observer, tag, info):
    if info is not None:
        info["timestamp"] = datetime.datetime.fromtimestamp(time.time()).strftime('[%H:%M:%S]')
        observer.notify("Scheduler", tag, info)
    else:
        observer.notify("Scheduler", tag, "")

def build_scheduler(scheduling_options, observer):
    """
    Factory function for scheduler building. Uses parameters to build up to 3 types of scheduler:
        - "Process/Parallel": A scheduler based on 'multiprocessing'. Can execute tasks in parallel.
        - "MPI/Parallel": This one uses MPI to execute tasks in parallel. Only usable with an MPI driver.
        - "Serial": A serial scheduler. Executes its tasks sequentially.

    @param scheduling_options: Parameters chunk containing the details of the scheduler to be created (mainly
    scheduler_type and number_of_processes).
    @param observer: The associated observer (can be None)

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
                                             },
                          'scheduling_started':{
                                             'function':send_message_to_observer,
                                             'kwargs':{
                                                       'observer':observer,
                                                       'tag':'Scheduler Starts'
                                                       }
                                             }
                          }

    scheduler_type = scheduling_options['scheduler_type']
    if scheduler_type == "Process/Parallel":
        max_processes = scheduling_options['number_of_processes'] if 'number_of_processes' in scheduling_options else multiprocessing.cpu_count()
        return ProcessParallelScheduler(max_processes, external_functions)

    elif scheduler_type == "MPI/Parallel":
        from pyscheduler.MPIParallelScheduler import MPIParallelScheduler # to avoid unneeded call to mpi_init
        return MPIParallelScheduler(share_results_with_all_processes=True, functions = external_functions)

    elif scheduler_type == "Serial":
        return SerialScheduler(external_functions)

    else:
        print "[ERROR][ClusteringExplorator::__init__] Not supported scheduler_type ( %s )"%(scheduler_type)
        exit()