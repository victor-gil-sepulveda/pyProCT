'''
Created on 04/02/2013

@author: victor
'''
import optparse
from pyproclust.protocol.protocolParameters import ProtocolParameters
from pyproclust.protocol.protocolImplementation import Protocol
import threading
from pyproclust.gui.observer.observer import Observer

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
            self.data_source.data_change_event.wait()
            print self.data_source.data
            self.data_source.data_change_event.clear()

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog [--nogui] script', version='1.0')
    
    #parser.add_option('-k', action="store", type='int', dest = "k", help="Number of clusters to get.",  metavar = "1")
    
    options, args = parser.parse_args()
    
    if(len(args)==0):
        parser.error("You need to specify the script to be executed.")
    
    json_script = args[0]
    
    parameters = ProtocolParameters.get_params_from_json(open(json_script).read())
    observer = Observer()
    
    cmd_thread = CmdLinePrinter(observer)
    cmd_thread.start()
    
    try:
        Protocol(observer).run(parameters)
    except:
        cmd_thread.stop()
        raise

    cmd_thread.stop()
