'''
Created on 21/01/2013

@author: victor
'''

import SimpleHTTPServer
import SocketServer
import logging
import webbrowser
import json
import hashlib
import time
from pyproclust.tools.scriptTools import create_directory
import os
from pyproclust.protocol.protocolImplementation import Protocol

def convert(my_input):
    if isinstance(my_input, dict):
        return {convert(key): convert(value) for key, value in my_input.iteritems()}
    elif isinstance(my_input, list):
        return [convert(element) for element in my_input]
    elif isinstance(my_input, unicode):
        return my_input.encode('utf-8')
    else:
        return my_input

if __name__ == '__main__':
    
    class ServerHandler(SimpleHTTPServer.SimpleHTTPRequestHandler):
        """
        Very simple implementation of a Request handler which only accepts POST requests.
        """
        
        def get_handlers(self):
            return {
                    "/run": self.run_handler,
                    "/save_params": self.save_params_handler,
                    "/file_exists": self.file_exists_handler,
                    "/create_directory": self.create_directory
                    }
        
        def file_exists_handler(self, data):
            data = convert(json.loads(data))
            print data
            self.wfile.write(json.dumps({"exists":os.path.exists(data['location']),
                                         "isfile":os.path.isfile(data['location']),
                                         "isdir":os.path.isdir(data['location'])}))
        
        def create_directory(self,data):
            data = convert(json.loads(data))
            print data
            try:
                success = create_directory(data['location'], ensure_writability = True)
                self.wfile.write(json.dumps({"done":success}))
            except:
                self.wfile.write(json.dumps({"done":False}))
           
        def run_handler(self, data):
            parameters = convert(json.loads(data))
            protocol = Protocol()
            protocol.run(parameters)
            pass
        
        def save_params_handler(self, data):
            data = convert(json.loads(data))
            create_directory("scripts")
            my_hash = hashlib.sha1()
            my_hash.update(str(time.time()))
            path = "scripts/"+my_hash.hexdigest()[:10]+".ppc"
            script_handler = open(path,"w")
            script_handler.write(json.dumps(data, sort_keys=False, indent=4, separators=(',', ': ')))
            script_handler.close()
            self.wfile.write('{"file_url":"'+path+'"}')
        
        def do_POST(self):
            fp= self.rfile
            data = fp.read(int(self.headers['Content-Length']))
            handle = self.get_handlers()[self.path]
            print "PATH", self.path
            handle(data)

    Handler = ServerHandler
    
    PORT = 8000
    httpd = SocketServer.TCPServer(("127.0.0.1", PORT), Handler)
    webbrowser.open("http://127.0.0.1:8000", new=0, autoraise=True)
    print "serving at port", PORT
    httpd.serve_forever()