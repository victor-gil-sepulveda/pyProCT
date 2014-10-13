'''
Created on Sep 3, 2014

@author: victor
'''
import copy

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
        
        # lexicographical ordering
        if self.source["source"] > cmp_string:
            return 1
        elif self.source["source"] < cmp_string:
            return -1
        else:
            return 0
    
    def get_path(self):
        return self.source["source"]
    
    def clone(self):
        return DataSource(copy.deepcopy(self.source))

    def add_info(self, key, info):
        self.source[key] = info
    
    def get_info(self, key):
        return self.source[key]
    
    def has_info(self, key):
        return key in self.source