'''
Created on 29/05/2012

@author: victor
'''
import numpy as np
from multiprocessing import Queue
from multiprocessing.process import Process
from pyproclust.clustering.analysis.analysisRunner import AnalysisRunner
from pyproclust.clustering.analysis.parallelAnalysis import analysis_to_parallel_analysis
from pyproclust.protocol.scheduling.processParallelScheduler import ProcessParallelScheduler

class ParallelAnalysisRunner(AnalysisRunner):
    def __init__(self,max_number_of_processes):
        '''
        Constructor
        '''
        self.analysis = []
        self.process_manager = ProcessParallelScheduler(max_number_of_processes,sleep_time = 30)
        self.analysis_id =  0
        self.analysis_queues = {}
        self.clusterings = []
        self.current_analysis = 0
        
    def add_analysis(self,analysis):
        """
        Adds one analysis type to the analysis queue.
        """
        panalysis = analysis_to_parallel_analysis(analysis)
        self.analysis.append(panalysis)
        self.analysis_queues[panalysis] = Queue()

    def run_analysis_for(self,clusterization):
        """
        Runs all the analysis we have for this clusterization.
        """
        self.clusterings.append(clusterization)
        
        def run_all_analysis_for_a_cluster():
            for a in self.analysis:
                a.run(clusterization, self.analysis_queues[a], self.current_analysis)
        func_kwargs = {}
        p = Process(target=run_all_analysis_for_a_cluster, name="Evaluation for " + str(self.current_analysis), kwargs=func_kwargs)
        p.description = "Analysis of: "+clusterization.details
        self.process_manager.add_process(p, p.name, [])
        self.current_analysis += 1

    
    def queue_unpacking(self,analysis,analysis_queue):
        """
        It parses a parallel queue with the results for a given analysis. After unpacking
        the data, it generates a string with the formatted results, and an array with the
        numerical results which will be returned.
        """
        res = analysis.results_string
        queue = analysis_queue
        analysis_results = []
        
        while not queue.empty():
            analysis_results.append(queue.get())
        analysis_results.sort()
        
        numerical_results = []
        for result in analysis_results:
            try:
                res += "%.3f\t"%(result[1])
                numerical_results.append(float(result[1]))
            except:
                res +=  result[1]+"\t"
                numerical_results.append(0)
            
        return res, numerical_results
    
    def numerical_results_normalization(self,numerical_results):
        """
        It gets the numerical results from queue_unpacking and normalizes them in a way that 
        the bigger one is 1
        """
        # make first element a float. If we can't, the result is not numerical
        try:
            numerical_results[0] = float(numerical_results[0])
            np_num_results = np.array(numerical_results)
            my_max = np.max(np.abs(np_num_results))
            if my_max != 0:
                np_num_results /= float(my_max)
        except:
            np_num_results = np.array(numerical_results)
            
        return np_num_results
    
    def generate_report(self):
        """
        Generates the final string and a repack of the numerical results, this time by cluster
         and not by analysis.
        """
        self.process_manager.consume()
        final_string,final_numerical_results = self.serial_gen_report(self.analysis,self.analysis_queues)
        return final_string, self.serial_result_packing(final_numerical_results,self.clusterings)
    
    def serial_gen_report(self,all_analysis,analysis_queues): 
        """
        Builds the final string and gathers all the analysis results (organized by analysis type).
        """    
        final_numerical_results = {}
        final_string = ""
        for a in all_analysis:
            res_string, numerical_results = self.queue_unpacking(a,analysis_queues[a])
            final_numerical_results[a] = self.numerical_results_normalization(numerical_results)
            final_string += res_string+"\n"
        return final_string,final_numerical_results
    
    def serial_result_packing(self,final_numerical_results,clusterings):
        """
        It reorganizes the numerical results by clustering and not by analysis type (meaning 
        that every cluster will be able to access to its analysis results).
        """    
        # Pack it
        results_pack = []
        for i in range(len(clusterings)):
            analysis_list = {}
            for a in final_numerical_results:
                analysis_list[a.name] = final_numerical_results[a][i]
            results_pack.append((clusterings[i],analysis_list))
        return results_pack
    