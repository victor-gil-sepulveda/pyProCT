"""
Created on 19/09/2012

@author: victor
"""
import os
import json
import shutil
import pyproct.tools.scriptTools as scripts_common
from pyproct.driver.observer.observable import Observable
from pyproct.driver.parameters import ProtocolParameters

class WorkspaceHandler(Observable):

    def __init__(self, workspace_parameters, observer):
        super(WorkspaceHandler,self).__init__(observer)

        self.parameters = ProtocolParameters(workspace_parameters.get_value("parameters", default_value = ProtocolParameters({
                                                                                                       "overwrite":True,
                                                                                                       "clear_after_exec":["tmp"]
                                                                                                  })))

        self.data = {
                      "results": os.path.join(workspace_parameters["base"], workspace_parameters.get_value("results", default_value="results")),
                      "tmp" : os.path.join(workspace_parameters["base"], workspace_parameters.get_value("tmp", default_value= "tmp")),
                      "clusters" : os.path.join(workspace_parameters["base"], workspace_parameters.get_value("clusters", default_value= "clusters")),
                      "matrix" : os.path.join(workspace_parameters["base"], workspace_parameters.get_value("matrix", default_value= "matrix"))
        }

    def __getitem__(self,key):
        return self.data[key]

    def __str__(self):
        return json.dumps(self.data, sort_keys=False, indent=4, separators=(',', ': '))

    def create_directories(self):
        """
        Recreates the workspace structure. Removes the old location if necessary.
        """
        self.notify("MSG","Creating workspace...")
        if self.parameters.get_value("overwrite", default_value = True)  :
            self.clear_directories(self.data.keys())

        for folder_key in self.data:
            scripts_common.create_directory(self.data[folder_key])

    def clear_directories(self, directory_keys):
        """
        Removes the directories given as parameters.

        @param directory_keys: The keys of the 'data' object which defined paths will be erased.
        """
        for folder_key in directory_keys:
            folder_path = self.data[folder_key]
            if os.path.exists(folder_path):
                shutil.rmtree(folder_path)
                self.notify("MSG","Removing %s ..."%folder_path)
#                 print "Removing %s ..."%folder_path

    def __enter__(self):
        self.create_directories()
        return self

    def __exit__(self, exception_type, exception_val, trace):
        self.clear_directories(self.parameters["clear_after_exec"])
