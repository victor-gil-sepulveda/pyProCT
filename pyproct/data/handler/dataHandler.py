"""
Created on 2/9/2014

@author: victor
"""
import glob
import os

class DataSource(object):
    
    def __init__(self, source):
        if isinstance(source, basestring):
            self.source = {"source": source}
        else:
            self.source = source
    
    def __cmp__(self, other):
        if isinstance(other, basestring):
            cmp_string = other
        else:
            cmp_string = other.source["source"]
        
        if self.source > cmp_string:
            return 1
        elif self.source < cmp_string:
            return -1
        else:
            return 0

class ElementRange(object):
    
    def __init__(self, start, end):
        self.start = start
        self.end = end
    
    def inflate(self):
        return range(self.start, self.end+1)
    
    def __contains__(self, item):
        return item>=self.start and item<=self.end
    
    def __iter__(self):
        for i in range(self.start, self.end+1):
            yield i

class SourceIterator(object):
    """
    Performs source expansion and yields data sources.
    """
    def __init__(self, source_list):
        self.source_list = []
        
        for source in source_list:
            if isinstance(source, basestring):
                self.source_list.append(DataSource(source))
            else:
                tmp_sources = glob.glob(source)
                for tmp_source in tmp_sources:
                    _, ext = os.path.splitext(tmp_source)
                    if ext == ".lst": #pyProCT file list
                        # load the dict
                        # TODO
                        # Then extend
                        self.source_list.extend([DataSource(list_source) for list_source in list_sources])
                    else:
                        self.source_list.append(DataSource(tmp_source))
                        
    def next(self):
        if len(self.source_list) == 0:
            raise StopIteration
        else:
            return self.source_list.pop(0)
        
    def __iter__(self):
        return self

class DataHandler(object):
    
    def __init__(self, params):
        """
        """
        self.elements = {}
        self.file_type = params["data"]["type"]
        for source in self.sources:
            element_range = self.loader.load(source)
            self.add_elements_with_same_source(element_range, source)
    
    
    def add_elements(self, elements, source_str):
        """
        """
        source = DataSource(source_str)
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
        
        