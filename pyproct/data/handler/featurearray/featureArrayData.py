'''
Created on 11/2/2015

@author: victor
'''
from pyproct.data.handler.data import Data
import numpy

class FeatureArrayData(Data):

    def __init__(self, data):
        self.data = data
        self.labels = data.keys()
        self.num_elements = len(self.data[self.data.keys()[0]])
        
    def get_number_of_elements(self):
        """
        Returns how much elements we have in this data package.
        """
        return self.num_elements
    
    def get_element(self, element_id):
        """
        Returns the datum associated to the element id passed as argument.
        """
        elements = []
        for label in self.labels:
            elements.append(self.data[label][element_id]) 
        return elements
    
    def get_all_elements(self):
        """
        Returns a matrix representation of all the elements.
        """
        return numpy.array([self.data[label] for label in self.labels]).T
    
    def get_feature_element(self, feature_label, element_id):
        return self.data[feature_label][element_id]
    
    def get_all_feature_elements(self, feature_label):
        return self.data[feature_label]
    
    def get_labels(self):
        return self.labels()