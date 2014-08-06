"""
Created on 04/02/2013

@author: victor
"""
from pyproct.driver.observer.accumulativeObserver import AccumulativeObserver

class Observable(object):
    """
    Definition of an observable class
    """

    def __init__(self, observer = None):
        if not observer is None:
            self.observer = observer
        else:
            self.observer = AccumulativeObserver()

    def notify(self, action, message):
        self.observer.notify(self.__class__.__name__, action, message)