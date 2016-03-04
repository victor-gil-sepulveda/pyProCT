"""
Created on Mar 4, 2016

@author: victor
"""

class OnlinePWHandler(object):
    
    def __init__(self, distance_provider):
        """
        Class constructor.
        
        :param distance_provider: A distance provider (something that gives distances by
        indexing a tuple e.g. distance = d_prov[item1,item2].
        """
        self.distance_matrix = distance_provider

    def save_matrix(self, matrix_path):
        """
        """
        print "[WARNING OnlinePWHandler::save_statistics] This operation will be skipped as it could be too costly."

    def save_statistics(self, matrix_base_path):
        """
        """
        print "[WARNING OnlinePWHandler::save_statistics] This operation will be skipped as it could be too costly."