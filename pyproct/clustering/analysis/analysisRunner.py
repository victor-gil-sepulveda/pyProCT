"""
Created on 29/05/2012

@author: victor
"""
def run_all_analysis_for_a_clustering(clustering_id, clustering, analysis):
    """
    Is the function to be run in parallel.
    
    @param clustering_id: Is the id of the clustering we are working with.
    
    @param clustering: A Clustering instance.
    
    @param analysis: A list of all the analysis we want to perform.s
    
    @param observer: An observer to communicate messages.
    """
    analysis_results = {}
    for a in analysis:
        analysis_results[a.name] =  a.run(clustering)
    return (clustering_id, analysis_results)

class AnalysisRunner():
    
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
        func_kwargs = {
                       "clustering_id": clustering_id,
                       "clustering":clustering_info[clustering_id]["clustering"],
                       "analysis":self.analysis
                       }
        
        task_name = "Evaluation of %s"%(clustering_id)
        description = "Evaluation of %s clustering (%s) with parameters: %s"%(clustering_info[clustering_id]["type"],
                                                                                 clustering_id,
                                                                                 clustering_info[clustering_id]["parameters"])
        
        self.scheduler.add_task(task_name = task_name, 
                                description = description, 
                                target_function = run_all_analysis_for_a_clustering, 
                                function_kwargs = func_kwargs,
                                dependencies = {})
        self.current_analysis += 1

                    
    def evaluate(self):
        """
        Runs all analysis for all clusterings in the clustering_info structure and recovers the results of this
        analsis to attach them into the clustering_info structure.
        """
        self.run_analysis_for_all_clusterings(self.clustering_info)
        
        results = self.scheduler.run()
        
        for (clustering_id, analysis_results) in results:
            self.clustering_info[clustering_id]["evaluation"] = analysis_results
    
