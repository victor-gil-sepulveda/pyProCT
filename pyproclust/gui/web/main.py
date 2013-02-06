'''
Created on 21/01/2013

@author: victor
'''

import SimpleHTTPServer
import SocketServer
import webbrowser
import json
import hashlib
import time
import os
from pyproclust.tools.scriptTools import create_directory
from pyproclust.protocol.protocolImplementation import Protocol
from pyproclust.tools.commonTools import convert_to_utf8
from pyproclust.protocol.protocolParameters import ProtocolParameters

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
            data = convert_to_utf8(json.loads(data))
            print data
            self.wfile.write(json.dumps({"exists":os.path.exists(data['location']),
                                         "isfile":os.path.isfile(data['location']),
                                         "isdir":os.path.isdir(data['location'])}))
        
        def create_directory(self,data):
            data = convert_to_utf8(json.loads(data))
            print data
            try:
                success = create_directory(data['location'], ensure_writability = True)
                self.wfile.write(json.dumps({"done":success}))
            except:
                self.wfile.write(json.dumps({"done":False}))
           
        def run_handler(self, data):
            parameters = ProtocolParameters.get_params_from_json(data)
            print parameters
            protocol = Protocol()
            protocol.run(parameters)
            pass
        
        def save_params_handler(self, data):
            data = convert_to_utf8(json.loads(data))
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