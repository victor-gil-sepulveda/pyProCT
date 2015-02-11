"""
Created on 5/9/2014

@author: victor
"""
from pyproct.data.handler.elementRange import ElementRange
import pyproct.tools.commonTools as common

class DataLoader(object):
    """
    Almost an interface for the data loader. A data loader must implement at least
    the close and load_data_from_source methods.
    """
    LOADER_TYPE = "none"
    
    def __init__(self, data_params):
        self.loaded_data = []
        self.number_of_elements = 0
        self.data_params = data_params
    
    def load(self, data_source):
        """
        Returns an ElementRange object with the loaded elements. For instance, 
        if we load a file containing 4 arrays and is the first time we load it,
        we will generate an ElementRange containing elements from 0 to 3 
        (0,1,2 and 3) representing each one of the loaded arrays. 
        
        TODO: when closed it cannot load anything else.
        """
        data, number_of_loaded_elements = self.load_data_from_source(data_source)
        self.loaded_data.append(data)
        current_number_of_elements = self.number_of_elements + number_of_loaded_elements
        e_range = ElementRange(self.number_of_elements, current_number_of_elements-1) 
        self.number_of_elements += len(e_range)
        return e_range
    
    def load_data_from_source(self, data_source):
        """
        Must return the data object and the number of loaded elements.
        """
        print "[ERROR DataLoader::load_data_from_source] This method must be overriden. Exiting..."
        exit()
    
    def close(self):
        """
        Must return the result of merging all the previously loaded data pieces as 
        a Data Object.
        If the number of loaded elements is 0 the program must exit (we 
        cannot perform any useful analysis without data!)  
        """
        if self.number_of_elements == 0:
            common.print_and_flush("[ERROR DataLoader:close] No loaded data. Exiting...\n")
            exit()