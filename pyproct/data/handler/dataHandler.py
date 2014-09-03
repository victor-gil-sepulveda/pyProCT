"""
Created on 2/9/2014

@author: victor
"""
from pyproct.data.handler.dataSource import DataSource
from pyproct.data.handler.sourceGenerator import SourceGenerator
from pyproct.tools.plugins import PluginHandler

class DataHandler(object):
    
    def __init__(self, params):
        """
        """
        self.elements = {}
        self.data_type = params["data"]["type"]
        self.sources = SourceGenerator(params["data"]["files"])
        
        loader_class = self.get_loader(self.file_type)
        loader = loader_class(params["data"])
        for source in self.sources:
            element_range = loader.load(source)
            self.add_elements_with_same_source(element_range, source)
        self.data = loader.close()
        
    def add_elements(self, elements, source):
        """
        """
        if not source in self.elements:
            self.elements[source] = elements
        else:
            self.elements[source].extend(elements)
        
    def get_elements(self, source_str):
        """
        """
        return self.elements[DataSource(source_str)].inflate()
        
    def get_source_of_element(self, element):
        """
        """
        sources = self.elements.keys()
        for source in sources:
            if element in self.elements[source]:
                return source
        return None
    
    def get_loader(self, data_type):
        # Get all available loaders
        available_loaders = PluginHandler.get_classes('pyproct.data.handler', 
                                                          selection_keyword = "DataLoader", 
                                                          skip_list = ["test"],
                                                          plugin_name = "data_loader")
        
        loaders = filter(lambda x: x.LOADER_TYPE == data_type, available_loaders)
        
        if len(loaders) == 0:
            print "[ERROR][DataHandler::get_loader] There is not a registered data loader for %s data type."%(data_type)
            exit()
        else:
            return loaders[0]
        
        
        
        