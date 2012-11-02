'''
Created on 16/08/2012

@author: victor
'''
from pyproclust.protocol.processPool import ProcessPool
import time

class SerialProcess():
    def __init__(self, t , name , kwargs):
        self.function = t
        self.name = name
        self.func_params = kwargs
    
    def run(self):
        kwargs = self.func_params
        self.function(**kwargs)
    
class SerialProcessPool(ProcessPool):
    def __init__(self):
        ProcessPool.__init__(self,0,0)
        self.sleep_time = 30
    
    def consume(self):
        
        process_names = self.processes.keys()
        
        while len( self.processes) > 0:
            print "--------"
            
            # Run a process
            process = self.processes.pop(process_names.pop())
            self.run_next_process(process)
            
            # Listing
            self.list_processes()
            
            # wait
            time.sleep(self.sleep_time) # 1 sec
        
        print "All processes have finished."
    
    def run_next_process(self,process):
        """
        Runs next process without accounting for dependencies.
        """
        process.run()
        self.finished.append(process.name)
        
    def add_process (self,process,name,dependency = []):
        if not name in self.processes.keys():
            self.processes[name] = process
            self.dependencies[name] = dependency
            self.num_added_processes += 1
        else:
            print "[Error ProcessPool::add_process] process already exists:",name
            exit()
    
    def add_process_internally(self,process_name,description,target_function,function_kwargs,dependency = []):
        process = SerialProcess(t = target_function, name = process_name, kwargs=function_kwargs)
        process.description = description
        self.add_process(process, process_name, dependency)
        