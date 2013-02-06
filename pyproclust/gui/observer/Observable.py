'''
Created on 04/02/2013

@author: victor
'''

class Observable(object):
    '''
    Definition of an observable class
    '''

    def __init__(self, observer):
        self.observer = observer
    
    def notify(self, action, message):
        self.observer.notify(self.__class__.__name__, action, message)