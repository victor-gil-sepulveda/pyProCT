"""
Created on 11/2/2015

@author: victor
"""
import os
import numpy
from pyproct.data.handler.dataLoader import DataLoader
from pyproct.data.handler.featurearray.featureArrayData import FeatureArrayData

class FeatureArrayDataLoader(DataLoader):
    """
    Loads an array of features from disk (numpy or txt).
    Layout must be row wise. First row can be the data labels.
    File must be numpy loadable. 
    """
    LOADER_TYPE = "features::array"
    
    def __init__(self, data_params):
        super(FeatureArrayDataLoader, self).__init__(data_params)
        self.feature_labels = []
        
    def close(self):
        """
        Prepares the merged array
        """
        super(FeatureArrayDataLoader, self).close()
        self.feature_labels = self.loaded_data[0].keys()
        
        merged_data = {}
        for data in self.loaded_data:
            if len(set(self.feature_labels)-set(data.keys())) != 0:
                print "[ERROR][FeatureArrayDataLoader::close] data sets with different number of labels were loaded."
                exit()
            for label in self.feature_labels:
                if label in merged_data:
                    merged_data[label].extend(data[label].tolist())
                else:
                    merged_data[label] = data[label].tolist()
        
        return FeatureArrayData(merged_data)

    def load_data_from_source(self, source):
        """
        :return: Prody's structure object with the loaded ensemble
        """
        _, ext = os.path.splitext(source.get_path())
        
        data = None
        if ext == ".npy":
            data = numpy.load(source.get_path())
        elif ext == ".txt":
            # It can have labels
            data = numpy.loadtxt(source.get_path())
        else:
            print "[ERROR][FeatureArrayDataLoader::load_data_from_source] file type not supported."%(source.get_path(),ext)
            exit()
            
        data_dic = {}
        labels = []
        if source.has_info("labels"):
            labels = source.get_info("labels")
            if len(labels) != len(data.T):
                print len(labels), len(data.T)
                print "[ERROR][FeatureArrayDataLoader::load_data_from_source] the number of labels is no equal to the number of columns."
                exit()
        else:
            labels = ["__feature_%d"%i for i in range(len(data.T))]
        
        for i, label in enumerate(labels):
            data_dic[label] = data.T[i]
        
        if source.has_info("load_only"):
            load_only = set(range(len(labels))) - set(source.get_info("load_only"))
            for column_id in  load_only:
                del data_dic[labels[column_id]]
        
        return  data_dic, len(data)
