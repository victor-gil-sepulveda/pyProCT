"""
Created on 2/9/2014

@author: victor
"""
from pyproct.data.handler.dataSource import DataSource
from pyproct.data.handler.sourceGenerator import SourceGenerator
from pyproct.tools.plugins import PluginHandler
from pyproct.data.handler.elementRange import ElementRange

class DataHandler(object):
    """
    Loads the data. Handles element ids and keeps a mapping with loaded datum. 
    Element ids are always in the range [0, number_of_datum]. The id for an element
    depends on the file load ordering (the order in the "file" array in the parameters
    and the order inside a file. 
    """
    def __init__(self, params, source_generator_class = SourceGenerator):
        """
        """
        self.elements = {}
        self.data_type = params["type"]
        self.sources = source_generator_class(params["files"]).source_list
        
        loader_class = self.get_loader(self.data_type)
        self.loader = loader_class(params)

        for source in self.sources:
            e_range = self.loader.load(source)
            self.add_elements(e_range, source)
        
        # Retrieve the data object
        self.data = self.loader.close()
    
    def get_number_of_elements(self):
        return self.data.get_number_of_elements()
    
    def get_data(self):
        return self.data
    
    def add_elements(self, elements, source):
        """
        """
        if not source in self.elements:
            self.elements[source] = elements
        else:
            self.elements[source].extend(elements)
        
    def get_elements_with_source(self, source_str):
        """
        """
        return [i for i in self.elements[DataSource(source_str)]]
    
    def get_all_elements(self):
        return ElementRange(0, self.get_number_of_elements()-1)
        
    def get_source_of_element(self, element):
        """
        """
        sources = self.elements.keys()
        for source in sources:
            if element in self.elements[source]:
                return source
        return None
    
    def get_loader(self, data_type):
        """
        Chooses the best loader for the type of data we have.
        """
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
        
    def __getstate__(self):
        """
        Pickling protocol (pickle)
        """
        state = {
                 "elements":self.elements,
                 "data_type":self.data_type,
                 "sources":self.sources,
                 "data": self.data
                 }
        return state
    
    def __setstate__(self, state):
        """
        Pickling protocol (unpickle)
        """
        self.elements = state["elements"]
        self.data_type = state["data_type"]
        self.sources = state["sources"]
        self.data = state["data"]