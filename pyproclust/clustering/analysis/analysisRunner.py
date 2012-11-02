'''
Created on 20/03/2012

@author: victor
'''

class AnalysisRunner(object):
    """
    This is a handler for the different kind of analysis we can do to a clusterization.
    It will run them all for our current clusterization and then return the tabbed string
    representing all the results (perfect for pasting in a datasheet).
    """
    def __init__(self):
        self.analysis = []
    
    def add_analysis(self,analysis):
        """
        Adds one analysis type to the analysis queue.
        """
        self.analysis.append(analysis)
    
    def run_analysis_for(self,clustering):
        """
        Runs all the analysis we have for this clustering.
        """
        for a in self.analysis:
            a.run(clustering)
    
    def generate_report(self):
        """
        Checks that all the analysis have the same ammount of results and 
        generates the report string.
        """
        ## Check that all the analysis have the same number of results,
        ## if not, it's time to abort.
        num_results = []
        for a in self.analysis:
            num_results.append(a.number_of_results)
        for i in range(len(self.analysis)-1):
            if num_results[i] != num_results[i+1]:
                print "[ERROR AnalysisRunner::generate_report] There was one analysis without the same ammount of results than the others."
        
        final_string = ""
        for a in self.analysis:
            final_string += a.results_string+"\n"
        return final_string