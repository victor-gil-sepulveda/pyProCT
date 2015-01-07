"""
Created on 26/11/2014

@author: victor
"""
import os
import numpy
from prody.measure.measure import calcCenter, calcDistance

class AtomicDistancesPostAction(object):
    KEYWORD = "atomic_distances"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, data_handler, workspace_handler, matrix_handler, generatedFiles):
        file_name = postprocessing_parameters.get_value("file", default_value = "atomic_distances") + ".csv"
        
        file_path = os.path.join( workspace_handler["results"], file_name)
        
        self.process_distances(file_path, clustering, postprocessing_parameters, data_handler, workspace_handler, matrix_handler, generatedFiles)
        
    def process_distances(self, file_path, clustering, postprocessing_parameters, data_handler, workspace_handler, matrix_handler, generatedFiles):
        handler = open(file_path,"w")
        distances =  calculate_selection_distances(postprocessing_parameters["distances"], data_handler)
        write_distances_file(handler, distances, clustering, postprocessing_parameters["distances"])
        handler.close()
        generatedFiles.append({
                               "description":"Atom distances file",
                               "path":os.path.abspath(file_path),
                               "type":"csv"
        })

def write_distances_file(handler, distances, clustering, distance_descriptors):
    
    for cluster in clustering.clusters:
        handler.write("cluster id, distance id, avg, std, elements in cutoff\n")
        for distance_id in distances:
            cluster_distances = distances[distance_id][cluster.all_elements]
            if "cutoff" in distance_descriptors[distance_id]:
                cutoff = float(distance_descriptors[distance_id]["cutoff"])
                print cluster_distances[cluster_distances < cutoff]
                handler.write("%s, %s, %.3f, %.3f, %d \n"%(cluster.id,
                                                           distance_id, 
                                                           numpy.mean(cluster_distances), 
                                                           numpy.std(cluster_distances),
                                                           len(cluster_distances[cluster_distances < cutoff])))
            else:
                handler.write("%s, %s, %.3f, %.3f \n"%(cluster.id,
                                                   distance_id, 
                                                   numpy.mean(distances[distance_id][cluster.all_elements]), 
                                                   numpy.std(distances[distance_id][cluster.all_elements])))
    
    handler.write("\n\n")
    
    for cluster in clustering.clusters:
        for distance_id in distances:
            handler.write("%s, %s, "%(cluster.id, distance_id) )
            for element in cluster.all_elements:
                handler.write("%.3f, "%(distances[distance_id][element]))
            handler.write("\n")

def calculate_selection_distances(distance_selection_objects, data_handler):
    """
    "distances":  {
           "dist_id_1":   {
                           "from": STRING::SELECTION,
                           "to": STRING::SELECTION,
                           "cutoff": FLOAT
                        },
                        ...
    }
    if cutoff is defined, only distances lower than this cutoff will be stored
    """
    distances = {}
    ensemble = data_handler.get_data().get_all_elements()
    
    for distance_id in distance_selection_objects:

        # Get selections for all frames
        try:
            from_selection = ensemble.select(distance_selection_objects[distance_id]["from"]) 
            to_selection = ensemble.select(distance_selection_objects[distance_id]["to"])
        except KeyError:
            print '[Error calculate_selection_distances ] "from" and "to" properties are mandatory. Skipping %s calculation...'%distance_id
            break
        
        # Check selections were ok
        if from_selection is None:
            print '"[Error calculate_selection_distances ] Unproductive "from" selection. Skipping %s calculation...'%distance_id
            break
        
        if to_selection is None:
            print '"[Error calculate_selection_distances ] Unproductive "to" selection. Skipping %s calculation...'%distance_id
            break
        
        distances[distance_id] = calcDistance(calcCenter(from_selection.getCoordsets()), 
                                              calcCenter(to_selection.getCoordsets()))
    
    return distances
    