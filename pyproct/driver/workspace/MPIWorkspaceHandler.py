"""
Created on 16/07/2014

@author: victor
"""
from pyproct.driver.workspace.workspaceHandler import WorkspaceHandler

class MPIWorkspaceHandler(WorkspaceHandler):

    def __init__(self, rank, allowed_rank, workspace_parameters, observer):
        """
        Constructor
        """
        super(MPIWorkspaceHandler,self).__init__(workspace_parameters, observer)
        self.disk_modif_allowed = (rank == allowed_rank)

    def create_directories(self):
        if self.disk_modif_allowed:
            super(MPIWorkspaceHandler,self).create_directories()

    def clear_directories(self, directory_keys):
        if self.disk_modif_allowed:
            super(MPIWorkspaceHandler,self).clear_directories(directory_keys)
