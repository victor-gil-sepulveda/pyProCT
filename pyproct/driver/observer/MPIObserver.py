"""
Created on 04/02/2013

@author: victor
"""
from pyproct.driver.observer.observer import ObservableMessage
from mpi4py import MPI
class MPIObserver(object):
    """
    Observer for MPI. It's replicated n times (where n is the number of processes), no possible race conditions
    occur except a possible race condition with stdout. Because of this, the observer accumulates all messages to
    finally write them into a file.
    """
    def __init__(self):
        self.data = None
        self.flag = False
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.messages = []

    def notify(self, actor, action, message):
        self.data = ObservableMessage(actor, action, message)
        self.messages.append(str(self.data))
        self.flag = True

    def get_data(self):
        return self.data

    def wait(self):
        while not self.flag:
            pass

    def clear(self):
        self.flag = False

    def __del__(self):
        "When destroyed, create a file with the contents."
        open("%03d.out"%self.rank,"w").write("".join([msg+"\n" for msg in self.messages]))