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
    
    def get_number_of_elements(self):
        """
        Returns how much elements we have in this data package.
        """
        print "[ERROR Data::get_number_of_elements] This method must be overriden. Exiting..."
        exit()
    
    def get_element(self, element_id):
        """
        Returns the datum associated to the element id passed as argument.
        """
        print "[ERROR Data::get_element] This method must be overriden. Exiting..."
        exit()
        
    def get_elements(self, element_list):
        """
        Returns an array with all the requested elements. Can be overriden to improve the 
        retrieving speed or changing the packaging style.
        """
        return [self.get_element(i) for i in element_list]
    
    def get_all_elements(self):
        """
        Returns a representation of all the elements.
        """
        print "[ERROR Data::get_all_elements] This method must be overriden. Exiting..."
        exit()
    