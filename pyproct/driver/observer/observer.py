"""
Created on 04/02/2013

@author: victor
"""
import threading
import copy

class ObservableMessage(object):

    def __init__(self, actor, action, message):
        self.contents = {
                         "actor":actor,
                         "action":action,
                         "message":message
        }

    def __str__(self):
        return "["+str(self.contents["actor"])+"]["+str(self.contents["action"])+"] "+str(self.contents["message"])

class Observer(object):
    """
    Observer class for the GUI version, thread safe.
    """

    def __init__(self):
        """
        Builds a new observer object.
        """
        #inits the semaphore
        self.semaphore = threading.Lock()
        self.data = None
        self.data_change_event = threading.Event()

    def notify(self, actor, action, message):
        self.semaphore.acquire()
        self.data = ObservableMessage(actor, action, message)
        self.semaphore.release()
        self.data_change_event.set()

    def acquire(self):
        self.semaphore.acquire()

    def release(self):
        self.semaphore.release()

    def get_data(self):
        self.semaphore.acquire()
        data = self.data
        self.semaphore.release()
        return copy.deepcopy(data)

    def wait(self):
        self.data_change_event.wait()

    def clear(self):
        self.data_change_event.clear()
