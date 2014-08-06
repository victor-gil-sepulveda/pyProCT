"""
Created on 07/02/2013

@author: victor
"""

class MeanObserver(object):
    def __init__(self):
        pass
    
    def notify(self, actor, action, message):
        print "I was notified of something by %s, but I'm not going to do anything about it."%actor
