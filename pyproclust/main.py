'''
Created on 04/02/2013

@author: victor
'''
import optparse
import threading
from pyproclust.driver.parameters import ProtocolParameters
from pyproclust.driver.observer.observer import Observer
from pyproclust.driver.driver import Driver
from pyproclust.driver.observer.MPIObserver import MPIObserver

class CmdLinePrinter(threading.Thread):
    
    def __init__(self,data_source):
        super(CmdLinePrinter, self).__init__()
        self._stop = threading.Event()
        self.data_source = data_source
        
    def stop(self):
        self._stop.set()
        self.data_source.notify("Main","Stop","Finished")

    def stopped(self):
        return self._stop.isSet()
    
    def run(self):
        while not self.stopped():
            self.data_source.wait()
            print self.data_source.get_data()
            self.data_source.clear()

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog [--mpi] [--gui] script', version='1.0')
    
    parser.add_option('--mpi', action="store_true",  dest = "use_mpi", help="Add this flag if you want to use MPI-based scheduling.")
    
    options, args = parser.parse_args()
    
    if(len(args)==0):
        parser.error("You need to specify the script to be executed.")
    
    json_script = args[0]
    
    parameters = None
    try:
        parameters = ProtocolParameters.get_params_from_json(open(json_script).read())
    except ValueError, e:
        print "Malformed json script."
        print e.message
        exit()
        
    observer = None
    if options.use_mpi:
        observer = MPIObserver()
    else:
        observer = Observer()
    
    cmd_thread = CmdLinePrinter(observer)
    cmd_thread.start()
    
    try:
        if not options.use_mpi:
            Driver(observer).run(parameters)
        else:
            from pyproclust.driver.mpidriver import MPIDriver
            MPIDriver(observer).run(parameters)
    except Exception, msg:
        print msg
        cmd_thread.stop()
        raise

    cmd_thread.stop()
