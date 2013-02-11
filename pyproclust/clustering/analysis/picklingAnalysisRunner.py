'''
Created on 29/05/2012

@author: victor
'''
import os
import tempfile
import pickle

def run_all_analysis_for_a_clustering(clustering_id, clustering, tmp_file_handler_descriptor, analysis):
    """
    Is the function to be run in parallel.
    
    @param clustering_id: Is the id of the clustering we are working with.
    
    @param clustering: A Clustering instance.
    
    @param tmp_file_handler_descriptor: The file handler of the temporary file created to hold the results.
    
    @param analysis: A list of all the analysis we want to perform.s
    """
    tmp_file_handler = os.fdopen(tmp_file_handler_descriptor,'w')
    analysis_results = {}
    for a in analysis:
        analysis_results[a.name] =  a.run(clustering)
    result = (clustering_id, analysis_results)
    pickle.dump(result, tmp_file_handler)
    del analysis_results
    tmp_file_handler.close()

class PicklingAnalysisRunner():
    
    def __init__(self, scheduler, parameters, clustering_info, populator):
        """
        Constructor
        
        @param scheduler: A Scheduler instance of any type to control the execution of analysis.
        
        @param parameters: General script parameters.
        
        @param clustering_info: a 'clustering_info' structure (dictionary that for each 'clustering_id' has a structure
        with all the info of one clustering. 
        
        @param populator: An 'AnalysisPopulator'-like instance to add the needed analysis.
        
        """
        self.analysis = populator.get_analysis_list()
        self.scheduler = scheduler
        self.evaluation_data = []
        self.current_analysis = 0
        self.clustering_info = clustering_info
    
    def add_analysis(self, analysis):
        """
        Simply appends one analysis to the analysis queue (the analysis we want to do).
        
        @param analysis: An 'Analysis' instance.
        """
        self.analysis.append(analysis)
       
    def run_analysis_for_all_clusterings(self, clustering_info):
        """
        Adds one process to the scheduler which will run all the analysis in the analysis queue for each
        of the clusterings in the clustering_info structure. The analysis are not executed instantaneously 
        but when the calculating function is used.
        
        @param clustering_info: The 'clustering_info' structure with the clustering to be analyzed.
        """
        for clustering_id in clustering_info:
            self.run_analysis_for_a_clustering(clustering_id, clustering_info)
        
    def run_analysis_for_a_clustering(self, clustering_id, clustering_info):
        """
        Adds one process to the scheduler which will run all the analysis in the analysis queue for one 
        clustering. The analysis are not executed instantaneously but when the calculating function is used.
        
        @param clustering_id: The clustering_id of the clustering we are going to work with.
        
        @param clustering_info: The 'clustering_info' structure with the clustering to be analyzed.
        """
        tmp_file_handler_descriptor , path = tempfile.mkstemp() 
        
        # Remember the path where we'll store the results
        self.evaluation_data.append(path)

        func_kwargs = {
                       "clustering_id": clustering_id,
                       "clustering":clustering_info[clustering_id]["clustering"],
                       "tmp_file_handler_descriptor":tmp_file_handler_descriptor,
                       "analysis":self.analysis
                       }
        
        process_name = "Evaluation of %s"%(clustering_id)
        description = "Evaluation of %s clustering (%s) with parameters: %s"%(clustering_info[clustering_id]["type"],
                                                                                 clustering_id,
                                                                                 clustering_info[clustering_id]["parameters"])
        
        self.scheduler.add_process(process_name, description, run_all_analysis_for_a_clustering, func_kwargs)
        self.current_analysis += 1

    def recover_evaluation_data(self, clustering_info):
        """
        Unpickles all the pickled results and builds the structures needed to start the bes clustering selection pipeline.
        """
        for path in self.evaluation_data:
            tmp_file_handler = open(path,"r")
            (clustering_id, analysis_results) = pickle.load(tmp_file_handler)
            tmp_file_handler.close()
            os.remove(path)
            clustering_info[clustering_id]["evaluation"] = analysis_results
                    
    def evaluate(self):
        """
        Runs all analysis for all clusterings in the clustering_info structure and recovers the results of this
        analsis to attach them into the clustering_info structure.
        """
        self.run_analysis_for_all_clusterings(self.clustering_info)
        self.scheduler.consume()
        self.recover_evaluation_data(self.clustering_info)
    
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
            for analysis_name in all_results.keys():
                analysis_dic[analysis_name] = all_results[analysis_name][i]
            repacked_results.append((clustering, analysis_dic))
        
        return repacked_results


    def gen_results_string(self, analysis_name, analysis_array):
        """
        Generates the report string for one analysis.
        """
        line =  analysis_name + "\t"
        for i in range(len(analysis_array)):
            line += str(analysis_array[i])+"\t"
        return line
       
        
    