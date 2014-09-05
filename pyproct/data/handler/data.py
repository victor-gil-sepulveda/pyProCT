"""
Created on 5/9/2014

@author: victor
"""

class Data(object):
    '''
    Interface for a data container. The data object holds all the data loaded from the source 
    files. It must implement one mandatory method ("get_element") and all the methods needed by
    any specialized matrix calculator, clustering validation indices or postprocessing analysis. 
    '''
    def __init__(self):
        pass
    
    def get_element(self, element_id):
        """
        Returns the datum associated to the element id passed as argument.
        """
        print "[ERROR Data::get_element] This method must be overriden. Exiting..."
        exit()