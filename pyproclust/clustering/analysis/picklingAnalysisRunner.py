'''
Created on 29/05/2012

@author: victor
'''
import os
import tempfile
import pickle
from pyproclust.driver.handlers.timerHandler import TimerHandler

def run_all_analysis_for_a_clustering(clustering_id, clustering, tmp_file_handler_descriptor, analysis):
    """
    Is the function to be run in parallel.
    
    @param clustering_id: Is the id of the clustering we are working with.
    
    @param clustering: A Clustering instance.
    
    @param tmp_file_handler_descriptor: The file handler of the temporary file created to hold the results.
    
    @param analysis: A list of all the analysis we want to perform.s
    
    @param observer: An observer to communicate messages.
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
    
    def __init__(self, scheduler, clustering_info, populator):
        """
        Constructor
        
        @param scheduler: A Scheduler instance of any type to control the execution of analysis.
        
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
    
