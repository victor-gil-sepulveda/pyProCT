"""
Created on 16/07/2014

@author: victor
"""

import inspect
import pkgutil
import importlib

class PostprocessingDriver(object):

    def __init__(self):
        self.available_actions = []

    def get_available_actions(self, root_package_p = 'pyproct.driver.sections.postprocessing.actions'):
        root_package = importlib.import_module(root_package_p)
        root_package_son_template_p = root_package_p+".%s"

        for pkg_info in pkgutil.walk_packages(root_package.__path__, onerror=lambda x: None):
            pckg_name = pkg_info[1]
            is_pkg = pkg_info[2]
            if is_pkg:
                self.get_available_actions(root_package_p="%s.%s"%(root_package_p, pckg_name))
            else:
                module_import_path = root_package_son_template_p%pckg_name
                module = importlib.import_module(module_import_path)
                for element_name, obj in inspect.getmembers(module):
                    if "PostAction" in element_name:
                        self.available_actions.append(obj())

    def run(self, clustering, postprocessing_parameters, trajectory_handler, workspace_handler, matrix_handler, generatedFiles):
        self.available_actions = []
        self.get_available_actions()
        for postprocessing_action in self.available_actions:
            postprocessing_action.run(clustering, postprocessing_parameters, trajectory_handler, workspace_handler, matrix_handler, generatedFiles)
