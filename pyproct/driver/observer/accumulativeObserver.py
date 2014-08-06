"""
Created on 07/02/2013

@author: victor
"""
class AccumulativeObserver(object):
    def __init__(self):
        self.messages = {}
    
    def notify(self, actor, action, message):
        if actor in self.messages.keys():
            self.messages[actor].append( {
                                    "action":action,
                                    "message":message
                                    })
        else:
            self.messages[actor] = [{
                                    "action":action,
                                    "message":message
                                    }]
        