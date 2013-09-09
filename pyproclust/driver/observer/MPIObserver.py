'''
Created on 04/02/2013

@author: victor
'''
from pyproclust.driver.observer.observer import ObservableMessage
class MPIObserver(object):
    '''
    Dummy observer.
    '''

    def __init__(self):
        self.data = None
        self.flag = False
    
    def notify(self, actor, action, message):
        self.data = ObservableMessage(actor, action, message)
        self.flag = True
    
    def get_data(self):
        return self.data
    
    def wait(self):
        while not self.flag:
            pass
    
    def clear(self):
        self.flag = False
    
