'''
Created on 27/11/2014

@author: victor
'''
import itertools
from prody.measure import calcDihedral
def sq_distance(a,b):
    return (a[0]-b[0])**2 + (a[1]-b[1]) + (a[2]-b[2])

def bonds_are_linked(b1,b2):
    return b1[0] == b2[0] or b1[1] == b2[1] or b1[0] == b2[1] or b1[1] == b2[0] 

def angles_share_bond(a1,a2):
    return a1[0] == a2[0] or a1[0] == a2[1] or a1[1] == a2[0] or a1[1] == a2[1]

def get_dihedral_indexes(dihedral):
    indexes = []
    for angle in dihedral:
        for bond in angle:
            indexes.append(bond[0], bond[1])
    return list(set(indexes))

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
        
        ref_selection = selection_coords[0]

        # Process bonds for reference frame (first)
        bonds = []
        sq_bond_distance = bond_distance**2
        for i in range(len(ref_selection)-1):
            for j in range(i+1, len(ref_selection)):
                if sq_distance(ref_selection[i],ref_selection[j]) <= sq_bond_distance:
                    bonds.append(tuple(sorted([i, j])))
        
        # Find angles
        angles = []
        for i in range(len(bonds)-1):
            for j in range(i+1, len(bonds)):
                if bonds_are_linked(bonds[i], bonds[j]):
                    angles.append(tuple(sorted([bonds[i], bonds[j]])))
        
        # Finally, find dihedrals
        dihedrals = []
        for i in range(len(bonds)-1):
            for j in range(i+1, len(bonds)):
                if angles_share_bond(angles[i], angles[j]):
                    dihedrals.append(tuple(sorted([angles[i], angles[j]])))
        
        # Now reorganize atoms in dihedrals so that 
        # they are consecutive and we can calculate the
        # actual dihedral angle
        r_dihedrals = []
        for dihedral in dihedrals:
            indexes = get_dihedral_indexes(dihedral)
            # Get permutation of minimum distance
            distances = []
            for perm in itertools.permutations(indexes):
                distances.append((sq_distance(ref_selection[perm[0]],ref_selection[perm[1]])+
                                  sq_distance(ref_selection[perm[1]],ref_selection[perm[2]])+
                                  sq_distance(ref_selection[perm[2]],ref_selection[perm[3]]),
                                  perm))
            # We will pick the one which summed distances is smaller 
            distances.sort()
            r_dihedrals.append(distances[0][1])
        
        all_angles = []
        for ref in selection_coords:
            #Calculate the angles for a ref
            angles = [] 
            for dihedral_indexes in r_dihedrals:
                atom1 = ref_selection[dihedral_indexes[0]]
                atom2 = ref_selection[dihedral_indexes[1]]
                atom3 = ref_selection[dihedral_indexes[2]]
                atom4 = ref_selection[dihedral_indexes[3]]
                
                angles.append(calcDihedral(atom1, atom2,  atom3, atom4))
            all_angles.append(angles)
            
        return all_angles
            
                
        
            