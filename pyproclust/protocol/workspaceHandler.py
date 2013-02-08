'''
Created on 19/09/2012

@author: victor
'''
import os
import json
import shutil
import pyproclust.tools.scriptTools as scripts_common
from pyproclust.protocol.observer.observable import Observable

class WorkspaceHandler(Observable):

    def __init__(self, parameters, observer):
        super(WorkspaceHandler,self).__init__(observer)
        
       
        self.data = {
                      "results": parameters["workspace"]["base"]+parameters["workspace"]["results"],
                      "tmp" : parameters["workspace"]["base"]+parameters["workspace"]["tmp"],
                      "clusterings" : parameters["workspace"]["base"]+parameters["workspace"]["clusterings"],
                      "matrix" : parameters["workspace"]["base"]+parameters["workspace"]["matrix"]
        }
        
        self.do_refinement = parameters["refinement"]["use"]
        
        if self.do_refinement:
            refinement_base = parameters["workspace"]["base"]+"/"+parameters["workspace"]["refinement_base"]
            self.data["refinement"]["pure_A"] = refinement_base+"/pure_A"
            self.data["refinement"]["pure_B"] = refinement_base+"/pure_B"
            self.data["refinement"]["mixed"] = refinement_base+"/mixed"
    
    def __getitem__(self,key):
        return self.data[key]
    
    def __str__(self):
        return json.dumps(self.data, sort_keys=False, indent=4, separators=(',', ': '))
    
    def create_directories(self, remove_existing = True):
        self.notify("","Creating workspace...")
        for path in self.data:
            if os.path.exists(self.data[path]) and remove_existing:
                self.notify("Removing Directory",self.data[path])
                shutil.rmtree(self.data[path])
            scripts_common.create_directory(self.data[path])
        
        if self.do_refinement:
            for path in self.data["refinement"]:
                if os.path.exists(self.data["refinement"][path]) and remove_existing:
                    self.notify("Removing Directory",self.data["refinement"][path])
                    shutil.rmtree(self.data["refinement"][path])
                scripts_common.create_directory(self.data["refinement"][path])
    
        
            
        