'''
Created on 21/01/2013

@author: victor
'''

import SimpleHTTPServer
import SocketServer
import logging
import cgi
import webbrowser
import json

if __name__ == '__main__':
    
    class ServerHandler(SimpleHTTPServer.SimpleHTTPRequestHandler):
        """
        Very simple implementation of a Request handler which only accepts POST requests.
        """
        
        def get_handlers(self):
            return {"/test": self.test_handler,
                    "/save_params": self.save_params_handler}
            
        def test_handler(self, data):
            print "DATA", data 
            logging.error(self.headers)
            self.wfile.write('{"caca":3}')
        
        def save_params_handler(self, data):
            print data
            print "DATA", json.loads(data,encoding='ASCII')
            logging.error(self.headers)
            self.wfile.write('{"caca":3}')
            
        def do_POST(self):
            fp= self.rfile
            data = fp.read(int(self.headers['Content-Length']))
            handle = self.get_handlers()[self.path]
            handle(data)

    Handler = ServerHandler
    
    PORT = 8000
    httpd = SocketServer.TCPServer(("127.0.0.1", PORT), Handler)
    webbrowser.open("http://127.0.0.1:8000", new=0, autoraise=True)
    print "serving at port", PORT
    httpd.serve_forever()