'''
Created on 31/01/2012

@author: victor
'''

class ImageStub(object):

    def __init__(self, w, h):
        self.matrix = []
        for r in range(h):
            row = []
            for i in range(w):
                row.append(0.0)
            self.matrix.append(row)
        del r
        del i
        
    def __getitem__(self, key):
        return self.matrix[key[0]][key[1]]
    
    def __setitem__(self, key,value):
        self.matrix[key[0]][key[1]] = value