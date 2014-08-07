"""
Created on 16/07/2014

@author: victor
"""

from pyproct.tools.plugins import PluginHandler

class PostprocessingDriver(object):
    
    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, trajectory_handler, workspace_handler, matrix_handler, generated_files):
        available_action_classes = PluginHandler.get_classes('pyproct.postprocess.actions', 
                                                                  "PostAction", 
                                                                  ["test","confSpaceComparison"],
                                                                  "postprocess")
        
        for postprocessing_action_class in available_action_classes:
            if postprocessing_action_class.KEYWORD in postprocessing_parameters:
                postprocessing_action_class().run(clustering, 
                                          postprocessing_parameters, 
                                          trajectory_handler, 
                                          workspace_handler, 
                                          matrix_handler, 
                                          generated_files)
