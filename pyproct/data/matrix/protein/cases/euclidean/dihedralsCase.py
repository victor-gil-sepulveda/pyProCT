'''
Created on 27/11/2014

@author: victor
'''
import itertools
import math
import numpy
import pyproct.tools.mathTools as mathTools

def bonds_are_linked(b1,b2):
    return b1[0] == b2[0] or b1[1] == b2[1] or b1[0] == b2[1] or b1[1] == b2[0] 

def angles_share_bond(a1,a2):
    return a1[0] == a2[0] or a1[0] == a2[1] or a1[1] == a2[0] or a1[1] == a2[1]

def get_dihedral_indices(dihedral):
    indices = []
    for angle in dihedral:
        for bond in angle:
            indices.append(bond[0])
            indices.append(bond[1])
    return list(set(indices))

def obtain_dihedral_angles(system_coords, bond_distance):
    """
    system_coords: coords for 1 frame
    """
    ref_selection = system_coords[0]
    
    # Process bonds for reference frame (first)
    bonds = []
    sq_bond_distance = bond_distance**2
    for i in range(len(ref_selection)-1):
        for j in range(i+1, len(ref_selection)):
            if mathTools.sq_distance(ref_selection[i], ref_selection[j]) <= sq_bond_distance:
                bonds.append(tuple(sorted([i, j])))
                
    print "DBG: Found %d bonds"%(len(bonds))
    
    # Find angles
    angles = []
    for i in range(len(bonds)-1):
        for j in range(i+1, len(bonds)):
            if bonds_are_linked(bonds[i], bonds[j]):
                angles.append(tuple(sorted([bonds[i], bonds[j]])))
    print "DBG: Found %d angles"%(len(angles))
    
    # Finally, find dihedrals
    dihedrals = []
    for i in range(len(angles)-1):
        for j in range(i+1, len(angles)):
            if angles_share_bond(angles[i], angles[j]):
                dihedrals.append(tuple(sorted([angles[i], angles[j]])))
    print "DBG: Found %d dihedrals"%(len(dihedrals))

    # Now reorganize atoms in dihedrals so that 
    # they are consecutive and we can calculate the
    # actual dihedral angle
    r_dihedrals = []
    for dihedral in dihedrals:
        indices = get_dihedral_indices(dihedral)
        # Get permutation of minimum distance
        distances = []
        for perm in itertools.permutations(indices):
            #print dihedral, perm
            distances.append(( mathTools.sq_distance(ref_selection[perm[0]],ref_selection[perm[1]])+
                               mathTools.sq_distance(ref_selection[perm[1]],ref_selection[perm[2]])+
                               mathTools.sq_distance(ref_selection[perm[2]],ref_selection[perm[3]]),
                              perm))
        # We will pick the one which summed distances is smaller 
        distances.sort()
        r_dihedrals.append(distances[0][1])
    
    all_angles = []
    for ref in system_coords:
        #Calculate the angles for a ref
        angles = [] 
        for dihedral_indexes in r_dihedrals:
            atom1 = ref[dihedral_indexes[0]]
            atom2 = ref[dihedral_indexes[1]]
            atom3 = ref[dihedral_indexes[2]]
            atom4 = ref[dihedral_indexes[3]]
            
            angles.append( mathTools.calc_dihedral(atom1, atom2,  atom3, atom4))
        all_angles.append(angles)
    return numpy.array(all_angles)

class DihedralEuclideanDistanceBuilder(object):

    def __init__(self, params):
        pass
    
    @classmethod
    def build(cls, data_handler, matrix_params):
        """
        If method = "euclidean_distance::ensemble" and "parameters::type" is "DIHEDRALS"
        We only need fit_selection and a "bond_distance" parameter.
        
        TODO: some var names are repeated and can be confusing
        """
        
        bond_distance = matrix_params.get_value("bond_distance", default_value=1.6)
        selection_coords = data_handler.get_data().getCalculationCoordinates()
        
        return obtain_dihedral_angles(selection_coords, bond_distance)
            
    @classmethod
    def calc_distances(cls, coordinates):
        distances = []
        for i in range(len(coordinates)-1): # number of frame
            for j in range(i+1, len(coordinates)): # frame j
                angular_increments = numpy.array([ mathTools.angular_distance(coordinates[i][k], coordinates[j][k]) for k in range(len(coordinates[i]))])
                # the distance is the sq root of the powd sum
                distances.append( math.sqrt((angular_increments**2).sum()))
        return numpy.array(distances)
        
            