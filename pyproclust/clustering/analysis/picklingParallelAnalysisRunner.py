'''
Created on 29/05/2012

@author: victor
'''
from pyproclust.protocol.processPool import ProcessPool
from multiprocessing.process import Process
from pyproclust.clustering.analysis.analysisRunner import AnalysisRunner
import tempfile
import pickle
from pyproclust.clustering.analysis.parallelAnalysisRunner import ParallelAnalysisRunner
import os

class PicklingParallelAnalysisRunner(ParallelAnalysisRunner,AnalysisRunner):
    
    def __init__(self,max_number_of_processes,sleep_time = 30):
        '''
        Constructor
        '''
        AnalysisRunner.__init__(self)
        self.process_manager = ProcessPool(max_number_of_processes,sleep_time)
        self.evaluation_data = []
        self.current_analysis = 0
     
    def add_analysis(self,analysis):
        """
        Simply appends one analysis type to be run.
        """
        AnalysisRunner.add_analysis(self,analysis) 
       
    def run_analysis_for(self,clustering):
        """
        Runs all the analysis we have for this clusterization.
        """
        tmp_file_handler_descriptor , path = tempfile.mkstemp() 
        self.evaluation_data.append(path)
        
        def run_all_analysis_for_a_cluster(clustering, tmp_file_handler_descriptor):
            tmp_file_handler = os.fdopen(tmp_file_handler_descriptor,'w')
            analysis_results = {}
            for a in self.analysis:
                analysis_results[a.name] =  a.run(clustering)
            result = (clustering,analysis_results)
            pickle.dump(result, tmp_file_handler)
            del analysis_results
            tmp_file_handler.close()
        
        func_kwargs = {"clustering":clustering,"tmp_file_handler_descriptor":tmp_file_handler_descriptor}
        process_name = "Evaluation nr. " + str(self.current_analysis)
        description = "Analysis of: " + clustering.details
        self.process_manager.add_process_internally(process_name, description, run_all_analysis_for_a_cluster, func_kwargs)
        
#        p = Process(target=run_all_analysis_for_a_cluster, name="Evaluation nr. " + str(self.current_analysis), kwargs=func_kwargs)
#        p.description = "Analysis of: " + clustering.details
#        self.process_manager.add_process(p, p.name, [])
        self.current_analysis += 1

    def gen_final_string_and_normalize(self, all_results):
        """
        As it name says, it will generate the final report string for all the analysis.
        It also normalizes the numerical arrays in the range [0,1].
        """
        final_string = ""
        for analysis_name in all_results:
            analysis_array = all_results[analysis_name] # Write the result string
            final_string += self.gen_results_string(analysis_name,analysis_array)+"\n" # Normalize all the numerical results
            # Hack!!!
            if not analysis_name in ("Number of clusters","Mean cluster size","Noise level"):
                all_results[analysis_name] = self.numerical_results_normalization(analysis_array)
            else:
                all_results[analysis_name] = analysis_array
        return final_string

    def repack_results(self, ordered_clusterings, all_results):
        """
        It repacks all the results in a way that each clustering has all its analysis results associated.
        This means that the result will be a tuple list where the first element is the clustering from which
        we did the analysis, and as second element it has a dictionary with elements of the type:
        "analysis name":"normalized value"
        Where the normalized value has been normalized over all the evaluations.
        """
        repacked_results = []
        for i in range(len(ordered_clusterings)):
            clustering = ordered_clusterings[i]
            analysis_dic = {}
            for ana_name in all_results.keys():
                analysis_dic[ana_name] = all_results[ana_name][i]
            repacked_results.append((clustering, analysis_dic))
        
        return repacked_results

    def recover_evaluation_data(self, ordered_clusterings, all_results):
        """
        Unpickles all the pickled results and builds structures to start the real evaluation.
        """
        for path in self.evaluation_data:
            tmp_file_handler = open(path,"r")
            result = pickle.load(tmp_file_handler)
            tmp_file_handler.close()
            os.remove(path)
            clustering, analysis_results = result
            ordered_clusterings.append(clustering)
            for analysis_result_name in analysis_results.keys():
                partial = analysis_results[analysis_result_name]
                try:
                    all_results[analysis_result_name].append(partial)
                except:
                    all_results[analysis_result_name] = [partial]

    def generate_report(self):
        """
        Generates the final string and a repack of the numerical all_results, this time by cluster
         and not by analysis.
        """
        self.process_manager.consume()
        
        ordered_clusterings = []
        all_results = {}
        self.recover_evaluation_data(ordered_clusterings, all_results)
        
        # Generate string and normalize
        final_string = self.gen_final_string_and_normalize(all_results)
        
        # Repack results
        repacked_results = self.repack_results(ordered_clusterings, all_results)    
        
        return final_string, repacked_results
    
    def gen_results_string(self, analysis_name, analysis_array):
        """
        Generates the report string for one analysis.
        """
        line =  analysis_name + "\t"
        for i in range(len(analysis_array)):
            line += str(analysis_array[i])+"\t"
        return line
       
        
    