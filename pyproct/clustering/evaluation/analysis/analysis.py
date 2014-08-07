"""
Created on 04/06/2012

@author: victor
"""
class Analysis(object):
    def __init__(self, name, analysis_function, other_params = None):
        """
        Creates one analysis object. It uses an 'analysis_function' which
        has at least a clustering as parameter and returns a string without any tab
        character.
        The analysis function can have another parameter. If we need more
        than one extra parameter it can be provided with an iterable.
        """
        self.name = name
        self.analysis_function = analysis_function # Function that will be called to do the calculation
        self.other_params = other_params

    def run(self, clustering):
        """
        Runs one analysis function and appends the result to the results string.
        It the result of the function is not an string, it converts the result to a 
        string.
        """
        if self.other_params:
            return self.analysis_function(clustering, self.other_params)
        else:
            return self.analysis_function(clustering)
    