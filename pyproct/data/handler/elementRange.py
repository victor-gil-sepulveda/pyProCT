"""
Created on Sep 3, 2014

@author: victor
"""

class ElementRange(object):
    
    def __init__(self, start, end):
        self.start = start
        self.end = end
    
    def __contains__(self, item):
        return item>=self.start and item<=self.end
    
    def __iter__(self):
        def xrange_generator(start, end):
            for i in xrange(self.start, self.end+1):
                yield i
        return xrange_generator(self.start, self.end)
            