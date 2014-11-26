"""
Created on 26/11/2014

@author: victor
"""
import os
from prody.measure.measure import calcCenter, calcDistance

class AtomicDistancesPostAction(object):
    KEYWORD = "atomic_distances"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, data_handler, workspaceHandler, matrixHandler, generatedFiles):
        file_name = parameters.get_value("file", default_value = "atomic_distances") + ".csv"
        handler = open(file_name,"w")
        distances =  calculate_selection_distances(postprocessing_parameters["distances"], data_handler)
        write_distances_file(handler, distances, clustering)
        handler.close()
        
        generatedFiles.append({
                               "description":"Atom distances file",
                               "path":os.path.abspath(file_name),
                               "type":"pdb"
        })

def write_distances_file(handler, distances):
    
    for cluster in clustering.clusters:
    for distance_id in distances:
        handler.write(distance_id+", ")
    
    for i in range(len)
        


def calculate_selection_distances(distance_selection_objects, data_handler):
    """
    "distances":  {
           "dist_id_1":   {
                           "from":"",
                           "to":""
                        },
                        ...
    }
    """
    distances = {}
    ensemble = data_handler.get_data().get_all_elements()
    
    for distance_id in distance_selection_objects:

        # Get selections for all frames
        try:
            from_selection = ensemble.select(distance_selection_objects[distance_id]["from"]) 
            to_selection = ensemble.select(distance_selection_objects[distance_id]["to"])
        except KeyError:
            print '[Error calculate_selection_distances ] "from" and "to" properties are mandatory. Skipping %s ...'%distance_id
            break
        
        # Check selections were ok
        if from_selection is None:
            print '"[Error calculate_selection_distances ] Unproductive "from" selection. Skipping %s ...'%distance_id
            break
        if to_selection is None:
            print '"[Error calculate_selection_distances ] Unproductive "to" selection. Skipping %s ...'%distance_id
            break
        
        distances[distance_id] = calcDistance(calcCenter(from_selection.getCoordsets()), 
                                              calcCenter(to_selection.getCoordsets()))
    
    return distances
    