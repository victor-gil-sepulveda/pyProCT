'''
Created on 06/06/2012

@author: victor
'''
from pyproclust.clustering.analysis.analysis import Analysis

def analysis_to_parallel_analysis(analysis):
    return ParallelAnalysis(analysis.name, analysis.analysis_function, analysis.other_params)

class ParallelAnalysis(Analysis):

    def __init__(self, name, analysis_function, other_params = None):
        Analysis.__init__(self,name, analysis_function, other_params)
        
    
    def run(self, clusterization, queue, clustering_id):
        """
        This version of the function is prepared to be used with 'multiprocessing' in a 
        concurrent execution. Results are then stored in a multiprocessing Queue object
        which is thread and process safe. With the clustering also a clustering_id must be
        provided (is different for each cluster) so that when the parallel processes end we are
        able to know which result goes with which cluster (executions of the analysis may not
        be in the sequential order).
        """
        if self.other_params:
            queue.put((clustering_id,self.analysis_function(clusterization,self.other_params)))
        else:
            queue.put((clustering_id,self.analysis_function(clusterization)))
        
        self.number_of_results += 1
