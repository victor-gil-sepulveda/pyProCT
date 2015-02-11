'''
Created on 1/12/2014

@author: victor
'''

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
from pyproct.data.matrix.protein.cases.euclidean.dihedralsCase import obtain_dihedral_angles
import numpy

class proteinDihedralPlotPostAction(object):
    KEYWORD = "protein_dihedral_2D_plot"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, data_handler, workspace_handler, matrix_handler, generatedFiles):
        file_name = postprocessing_parameters.get_value("file", default_value = "dihedral_plot") + ".png"
        file_path = os.path.join( workspace_handler["results"], file_name)
        
        if not "dihedrals" in postprocessing_parameters:
            raise KeyError("[Error proteinDihedralPlot::run] \"dihedrals\" property is mandatory.")
        
        dihedral_descriptor = postprocessing_parameters.get_value("dihedrals")
        
        x_selection_s = dihedral_descriptor["x"]
        y_selection_s = dihedral_descriptor["y"]
        
        x_selection = data_handler.get_data().getSelectionCoordinates(x_selection_s)
        y_selection = data_handler.get_data().getSelectionCoordinates(y_selection_s)
        
        if x_selection is None or y_selection is None:
            raise KeyError("[Error proteinDihedralPlot::run] Unproductive selections for dihedral calculation.")
        
        x_values = obtain_dihedral_angles(x_selection, 1.6)
        y_values = obtain_dihedral_angles(y_selection, 1.6)
        
        # check values
        if len(x_values[0]) > 1 or len(y_values[0]) > 1:
            print x_values, y_values
            raise KeyError("[Error proteinDihedralPlot::run] Selection defined more than one dihedral.")
        
        x_values = x_values.flatten()
        y_values = y_values.flatten()
        
        colors = iter(cm.rainbow(numpy.linspace(0, 1, len(clustering.clusters))))
        for cluster in clustering.clusters:
                plt.scatter(x_values[cluster.all_elements], y_values[cluster.all_elements],color=next(colors))
        
        plt.savefig(file_path)
        plt.close()
        
        generatedFiles.append({
                               "description":"Dihedrals plot",
                               "path":os.path.abspath(file_path),
                               "type":"image"
        })
