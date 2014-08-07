"""
Created on 04/02/2013

@author: victor
"""
import time
import functools

class Timer(object):
    def __init__(self):
        self.start_t = 0.0
        self.end_t = 0.0
        self.stopped = True

    def start(self):
        self.stopped = False
        self.start_t = time.time()
        return self

    def stop(self):
        self.stopped = True
        self.end_t = time.time()
        return self

    def get_elapsed_time(self):
        if self.stopped :
            return self.end_t - self.start_t
        else:
            return time.time() - self.start_t

class TimerHandler(object):

    def __init__(self):
        self.named_timers = {}
        self.timer_list = []

    def start(self,name):
        self.named_timers[name] = Timer()
        self.named_timers[name].start()
        self.timer_list.append((name,self.named_timers[name]))

    def stop(self, name):
        self.named_timers[name].stop()

    def get_elapsed(self):
        timers = []
        for name, timer in self.timer_list:
            timers.append({
                          "name":name,
                          "elapsed":timer.get_elapsed_time()
            })
        return timers

    def __str__(self):
        rep = ""
        for timer in self.get_elapsed():
            rep += "%s : %.3f\n"%(timer["name"], timer["elapsed"])
        return rep

    def __repr__(self):
        return self.__str__()

class timed_method(object):
    def __init__(self, timer, alias):
        self.timer = timer
        self.alias = alias

    def __call__(self, f):
        def timer_wrap(*args, **kwargs):
            self.timer.start(self.alias)
            result = f(*args, **kwargs)
            self.timer.stop(self.alias)
            return result
        functools.update_wrapper(timer_wrap, f)
        return timer_wrap
