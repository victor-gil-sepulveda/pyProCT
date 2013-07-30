'''
Created on 01/07/2013

@author: victor
'''
class MPIParallelScheduler(object):
    
    def __init__(self, number_of_processes, sleep_time = 30, observer = None):
    
    def consume(self):
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        all_processes = comm.bcast(self.processes, root=0)
        while not len(self.finished) == self.num_added_processes:
            pass
        MPI.
        
    def list_processes(self):
        pass
    
    def run_next_process(self):
        pass
    
    def get_not_dependant_process(self):
        pass
    
    def add_process(self, process_name, description, target_function, function_kwargs, dependencies = []):
        if not process_name in self.processes.keys():
            process = Process(target = target_function, name = process_name, kwargs=function_kwargs)
            process.description = description
            self.processes[process_name] = process
            self.dependencies[process_name] = dependencies
            self.num_added_processes += 1
        else:
            print "[Error ProcessPool::add_process] process already exists:", process_name
            exit()
        
    def next_process_id (self):
        return self.num_added_processes