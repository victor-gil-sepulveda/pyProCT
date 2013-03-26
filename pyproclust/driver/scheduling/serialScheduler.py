'''
Created on 16/08/2012

@author: victor
'''
import time
from pyproclust.protocol.observer.observable import Observable
from pyproclust.protocol.scheduling.processParallelScheduler import ProcessParallelScheduler

class SerialProcess(Observable):
    def __init__(self, t, name, kwargs, observer = None):
        Observable.__init__(self, observer)
        self.function = t
        self.name = name
        self.func_params = kwargs
    
    def run(self):
        kwargs = self.func_params
        self.notify("[Process Begin]", {
                                        "process": self.name,
        })
        self.function(**kwargs)
        self.notify("[Process End]", {
                                        "process": self.name,
        })
    
class SerialScheduler(ProcessParallelScheduler):
    def __init__(self, sleep_time = 5, observer = None):
        ProcessParallelScheduler.__init__(self,0,0, observer)
        self.sleep_time = sleep_time
    
    def consume(self):
        process_names = sorted(self.processes.keys())
        while len( self.processes) > 0:
            # Run a process
            self.run_next_process(process_names)
            
            # Listing
            self.list_processes()
            
            # wait
            time.sleep(self.sleep_time) # 1 sec
        
        self.notify("[Ended]","")
    
    def run_next_process(self, process_names):
        """
        Runs next process without accounting for dependencies.
        """
        process = self.processes.pop(process_names.pop())
        process.run()
        self.finished.append(process.name)
    
    def add_process(self, process_name, description, target_function, function_kwargs):
        if not process_name in self.processes.keys():
            process = SerialProcess(t = target_function, name = process_name, kwargs=function_kwargs, observer = self.observer)
            process.description = description
            self.processes[process_name] = process
            self.num_added_processes += 1
        else:
            print "[Error ProcessPool::add_process] process already exists:", process_name
            exit()
        