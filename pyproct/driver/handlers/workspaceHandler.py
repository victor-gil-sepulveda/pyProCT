'''
Created on 19/09/2012

@author: victor
'''
import os
import json
import shutil
import pyproct.tools.scriptTools as scripts_common
from pyproct.driver.observer.observable import Observable

class WorkspaceHandler(Observable):

    def __init__(self, workspace_parameters, observer):
        super(WorkspaceHandler,self).__init__(observer)

        self.data = {
                      "results": os.path.join(workspace_parameters["base"], workspace_parameters["results"] if "results" in workspace_parameters else "results"),
                      "tmp" : os.path.join(workspace_parameters["base"], workspace_parameters["tmp"] if "results" in workspace_parameters else "tmp"),
                      "clusters" : os.path.join(workspace_parameters["base"], workspace_parameters["clusters"] if "clusters" in workspace_parameters else "clusters"),
                      "matrix" : os.path.join(workspace_parameters["base"], workspace_parameters["matrix"] if "matrix" in workspace_parameters else "matrix")
        }

    def __getitem__(self,key):
        return self.data[key]

    def __str__(self):
        return json.dumps(self.data, sort_keys=False, indent=4, separators=(',', ': '))

    def create_directories(self, remove_existing = True):
        self.notify("MSG","Creating workspace...")
        for path in self.data:
            if os.path.exists(self.data[path]) and remove_existing:
                self.notify("Removing Directory",self.data[path])
                shutil.rmtree(self.data[path])
            scripts_common.create_directory(self.data[path])



