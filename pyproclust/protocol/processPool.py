'''
Created on 16/05/2012

@author: victor
'''
from multiprocessing import Process
import time
from threading import Thread, Lock

class ProcessLauncher(Thread):
    def __init__(self,process,update_guard, processes_dict, running_processes_list, finished_processes_list):
        Thread.__init__(self)
        self.process = process
        self.update_guard = update_guard
        self.processes_dict_proxy = processes_dict
        self.currently_running_processes_list_proxy = running_processes_list
        self.finished_processes_proxy = finished_processes_list
        
    def run(self):
        # Update lists (process will start to run, and will no longer be available
        # for execution (so we delete it from the processes dictionary)
        self.update_guard.acquire()
        self.currently_running_processes_list_proxy.append(self.process)
        del self.processes_dict_proxy[self.process.name]
        self.update_guard.release()
        
        # Start the process
        self.process.start()
        # Wait for termination
        self.process.join()
        
        #Update lists
        self.update_guard.acquire()
        self.currently_running_processes_list_proxy.remove(self.process)
        self.finished_processes_proxy.append(self.process.name)
        self.update_guard.release()
    
class ProcessPool(object):  
    
    def __init__(self,max_processes_at_the_same_time, sleep_time = 30):
        self.finished = []
        self.processes = {}
        self.dependencies = {}
        self.currently_running_processes = []
        self.max_processes = max_processes_at_the_same_time
        self.num_added_processes = 0
        self.sleep_time = sleep_time
        self.thread_guard = Lock()
    
    def consume(self):
        while not len(self.finished) == self.num_added_processes:
            print "--------"
            
            # Run a process
            self.run_next_process()
            
            # Listing
            self.list_processes()
            
            # wait
            time.sleep(self.sleep_time) # 1 sec
        
        print "All processes have finished."
    
    def list_processes(self):
        self.thread_guard.acquire()
        for p in self.currently_running_processes:
            print "Running:",p.name
        for p in self.processes:
            print "Idle: ", self.processes[p].name
        for p in self.finished:
            print "Ended: ", p
        self.thread_guard.release()
    
    def run_next_process(self):
        self.thread_guard.acquire()
        process = self.get_not_dependant_process()
        self.thread_guard.release()
        
        if process == None:
            self.thread_guard.acquire()
            if len(self.finished) == len(self.processes.keys()):
                print "No more processes to launch."
            self.thread_guard.release()
        else:
            self.thread_guard.acquire()
            number_of_processes_being_executed = len(self.currently_running_processes)
            self.thread_guard.release()
            
            if number_of_processes_being_executed < self.max_processes:
                print "Launching process:", process.description,".",self.max_processes-len(self.currently_running_processes)-1,"processes left."
                ProcessLauncher(process, self.thread_guard, self.processes, self.currently_running_processes,  self.finished).start()
            
    def get_not_dependant_process(self):
        for process_name in self.processes.keys():
            if len(self.dependencies[process_name]) == 0:
                return self.processes[process_name]
            else:
#                print "No independent process."
                found_dependency  = False
                for d in self.dependencies[process_name]:
                    if not d in self.finished :
                        found_dependency = True
                        print "Dependendy found (",process_name ,"):",d
                if not found_dependency:
                    return self.processes[process_name]
#        print "No suitable process."
        return None
    
    def add_process (self,process,name,dependency = []):
        if not name in self.processes.keys():
            self.processes[name] = process
            self.dependencies[name] = dependency
            self.num_added_processes += 1
        else:
            print "[Error ProcessPool::add_process] process already exists:",name
            exit()
            
    def add_process_internally(self,process_name,description,target_function,function_kwargs,dependency = []):
        process = Process(target = target_function, name = process_name, kwargs=function_kwargs)
        process.description = description
        self.add_process(process, process_name, dependency)
        
    def next_process_id (self):
        return self.num_added_processes
    
# Some complex Tests
if __name__ == "__main__":
    def test_without_params():
        print "Sleeping 6 seconds"
        time.sleep(6)
        print "Done" 
    
    def test_with_params(seconds):
        print "Sleeping %d seconds"%(seconds)
        time.sleep(seconds)
        print "Done"
    
    p1 = Process(target = test_without_params,name = "Process 1")
    p1.description = "First process"
    p2 = Process(target = test_with_params,name = "Process 2",kwargs={"seconds":15})
    p2.description = "Second process"
    p3 = Process(target = test_with_params,name = "Process 3",kwargs={"seconds":10})
    p3.description = "Third process"

    # Without dependencies!
#    pool = ProcessPool(2)
#    pool.add_process(p1, p1.name)
#    pool.add_process(p2, p2.name)
#    pool.add_process(p3, p3.name)
#    pool.consume()
    
    
    # And with dependencies!
    pool = ProcessPool(2)
    pool.add_process(p2, p2.name,[p1.name])
    pool.add_process(p1, p1.name)
    pool.add_process(p3, p3.name,[p1.name])
    # This may obligue p1 to be executed and finished alone
    pool.consume()