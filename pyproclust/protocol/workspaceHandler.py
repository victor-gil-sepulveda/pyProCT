'''
Created on 19/09/2012

@author: victor
'''
import pyproclust.tools.scriptTools as scripts_common
import pyproclust.tools.commonTools as common

class WorkspaceHandler(object):

    def __init__(self,protocol_parameters):
        
        self.results_path = protocol_parameters.getWorkspacePathFor("results")
        self.matrix_path = protocol_parameters.getWorkspacePathFor("matrix analysis")
        self.tmp_path = protocol_parameters.getWorkspacePathFor("tmp")
        self.clusterings_path = protocol_parameters.getWorkspacePathFor("clusterings")
        self.do_refinement = protocol_parameters.do_refinement
        if self.do_refinement:
            self.__refinement_base = protocol_parameters.getWorkspacePathFor("refinement_base")
            self.refinement_pure_A = self.__refinement_base+\
                            protocol_parameters.getWorkspacePathFor("pure_A",add_working_dir = False)
            self.refinement_pure_B = self.__refinement_base+\
                            protocol_parameters.getWorkspacePathFor("pure_B",add_working_dir = False)
            self.refinement_mixed = self.__refinement_base+\
                            protocol_parameters.getWorkspacePathFor("mixed",add_working_dir = False)
         
    def create_directories(self):
        common.print_and_flush( "Creating workspace...")
        scripts_common.create_directory(self.results_path)
        scripts_common.create_directory(self.matrix_path)
        scripts_common.create_directory(self.tmp_path)
        scripts_common.create_directory(self.clusterings_path)
        if self.do_refinement:
            scripts_common.create_directory(self.refinement_pure_A)
            scripts_common.create_directory(self.refinement_pure_B)
            scripts_common.create_directory(self.refinement_mixed)
        common.print_and_flush( "Done\n")
        