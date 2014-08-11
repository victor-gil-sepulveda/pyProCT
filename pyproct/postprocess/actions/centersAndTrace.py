"""
Created on 16/07/2014

@author: victor
"""
import json
import os
import numpy
import matplotlib.cm as cm

class CentersAndTracePostAction(object):
    KEYWORD = "centers_and_trace"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):
        centers_path, centers_contents = generate_selection_centers_file(clustering,
                                                                        workspaceHandler,
                                                                        trajectoryHandler)

        open(centers_path,"w").write(json.dumps(centers_contents,
                                  sort_keys=False,
                                  indent=4,
                                  separators=(',', ': ')))
        generatedFiles.append({
                                    "description":"Centers of the selection used to calculate distances",
                                    "path":os.path.abspath(centers_path),
                                    "type":"text"
        })

                
def calculate_bounding_box(coordinates, backbone_trace = []):
    """
    Calculates the bounding box that encloses all the atoms of two arrays.
    @param coordinates: An array containing coordinates for a set of conformations (3D array where dim 0 is the
    conformation, dim 1 is the atom and dim 2 the coordinate)
    @param backbone_trace: If used (and if it has length > 0) is a 2D array where dim 0 is the atom and dim 1 the
    coordinate.
    @return: A 3-tuple containing:
    - The bounding box corners.
    - The center of the b.b. .
    - The corner with maximum coordinates values.
    """
    coords = numpy.array(coordinates)
    if(len(backbone_trace)>0):
        [max_x,max_y,max_z] = numpy.max([numpy.max(numpy.max(coords,1),0).tolist()]+[numpy.max(backbone_trace,0).tolist()],0)
        [min_x,min_y,min_z] = numpy.min([numpy.min(numpy.min(coords,1),0).tolist()]+[numpy.min(backbone_trace,0).tolist()],0)
    else:
        [max_x,max_y,max_z] = numpy.max(numpy.max(coords,1),0)
        [min_x,min_y,min_z] = numpy.min(numpy.min(coords,1),0)

    center = numpy.array([min_x,min_y,min_z]) + ((numpy.array([max_x,max_y,max_z])-numpy.array([min_x,min_y,min_z])) /2.)
    return ([[max_x, max_y, max_z],
            [max_x, max_y, min_z],
            [max_x, min_y, max_z],
            [max_x, min_y, min_z],
            [min_x, max_y, max_z],
            [min_x, max_y, min_z],
            [min_x, min_y, max_z],
            [min_x, min_y, min_z]], center.tolist(), [max_x,max_y,max_z])

def generate_CA_or_P_trace(trajectoryHandler, backbone_atoms_selection = "name CA P"):
    """
    Gets the coordinates of the atoms forming the backbone of a protein. By default We consider the CA atoms in
    proteins and P atoms in DNA/RNA, but of course is an arbitrary choice.
    @param trajectoryHandler: Is the project's trajectory handler.
    @param backbone_atoms_selection: Selection describing the atoms that form part of the trace.
    @return: A list containing the atom positions of the backbone trace (ordered).
    """
    coordsets = numpy.array([])
    try:
        # Only get first frame of the selection
        coordsets = trajectoryHandler.getMergedStructure().select(backbone_atoms_selection).getCoordsets()[0]
    except:
        print "[ERROR visualizationTools::generate_CA_or_P_trace] Impossible to get coordinates for trace"
    return coordsets.tolist()

def generate_selection_centers_file(clustering, workspaceHandler, trajectoryHandler):
    # TODO: Superpose and center coords (or getting already superposed confs)
    #########################

    #########################
    centers_path = os.path.join(workspaceHandler["results"], "selection_centers.json")
    ligand_coords = trajectoryHandler.getCalculationCoordinates()

    centers_contents={}
    centers = []

    # Calculate trace
    centers_contents["backbone_trace"] = generate_CA_or_P_trace(trajectoryHandler)

    # Get Bounding Box
    (centers_contents["bounding_box"] ,
    centers_contents["bounding_box_center"],
    centers_contents["bounding_box_corner"]) = calculate_bounding_box( ligand_coords.tolist() ,
                                                                      centers_contents["backbone_trace"])

    # Colors iterator
    colors = iter(cm.rainbow(numpy.linspace(0, 1, len(clustering.clusters))))

    # calculate per cluster centers for selection (and prototype center)
    centers_contents["points"] = {}
    for cluster in clustering.clusters:
        centers = []
        for element in cluster.all_elements:
            coords = ligand_coords[element]
            centers.append(list(coords.mean(0)))
        centers_contents["points"][cluster.id] = {}
        centers_contents["points"][cluster.id]["prototype"] = list(ligand_coords[cluster.prototype].mean(0))
        centers_contents["points"][cluster.id]["centers"] = centers
        centers_contents["points"][cluster.id]["color"] = list(next(colors))[0:3]

    centers_contents["percents"] = {}
    for cluster in clustering.clusters:
        centers_contents["percents"][cluster.id] = "%.3f"%((len(cluster.all_elements) / float(clustering.total_number_of_elements))*100)

    return centers_path, centers_contents
