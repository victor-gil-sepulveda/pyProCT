'''
Created on 19/09/2012

@author: victor
'''
import pyproclust.tools.scriptTools as scripts_common
import pyproclust.tools.commonTools as common
import os
import json

class WorkspaceHandler(object):

    def __init__(self, parameters):
        self.data = {
                      "results": parameters["workspace"]["base"]+parameters["workspace"]["results"],
                      "tmp" : parameters["workspace"]["base"]+parameters["workspace"]["tmp"],
                      "clusterings" : parameters["workspace"]["base"]+parameters["workspace"]["clusterings"]
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
        common.print_and_flush( "Creating workspace...\n")
        for path in self.data:
            if os.path.exists(self.data[path]) and remove_existing:
                common.print_and_flush(self.data[path]+" exists, removing.\n")
                os.rmdir(self.data[path])
            scripts_common.create_directory(self.data[path])
        
        if self.do_refinement:
            for path in self.data["refinement"]:
                if os.path.exists(self.data["refinement"][path]) and remove_existing:
                    common.print_and_flush(self.data["refinement"][path]+" exists, removing.\n")
                    os.rmdir(self.data["refinement"][path])
                scripts_common.create_directory(self.data["refinement"][path])
        
        common.print_and_flush( "Done\n")
    
        
            
        