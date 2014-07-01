'''
Created on 27/06/2014

@author: victor
'''
import numpy
import itertools
import pyRMSD.RMSDCalculator
import pyRMSD.condensedMatrix
from pyproct.tools.prodyTools import removeAllCoordsetsFromStructure

class ChainMappingRMSDMatrixCalculator:
    """
    Some structures can be formed of repeated units that can be arbitrarily place. This class tries to reorder this units (chains)
    in order to calculate the minimum possible RMSD. Ex. Given 3 structures formed of three equal domains:


                        B
        Structure 1:  A   C

                        A
        Structure 2:  B   C

                        B
        Structure 3:  A   C

    With resulting coordsets (using first ordering as reference):

             1           2          3
      1  [A][B][C]   [B][C][A]   [B][A][C]

      2              [A][B][C]   [C][B][A]

      3                          [A][B][C]

    This class calculates the RMSD with the mapping that produces the lower RMSD. Ex.
    for 1 -> 2 the best mapping is A -> B, B -> C and C -> A
    for 1 -> 3 the best mapping is A -> B, B -> A and C -> C
    for 2 -> 3 the best mapping would be A -> C, B -> B and C -> A

    Currently only works if all chains have exactly the same length. If chains are of different lengths
    then the result is not well defined (padding would be needed).


    Improvement: Only do permutations of chains with the same length! -> Trivial solution with padding
    """
    def __init__(self):
        pass

    @classmethod
    def getStructureChains(cls, pdb_structure, in_chain_selection):
        """

        """
        chains = pdb_structure.select(in_chain_selection).copy().getHierView()
        number_of_models = pdb_structure.numCoordsets()
        structure_chains = []
        for i in range(number_of_models):
            struct_chain = {}
            for chain in chains:
                struct_chain[chain.getChid()] = chain.getCoordsets(i)
            structure_chains.append(struct_chain)
        return structure_chains

    @classmethod
    def reorderCoordinates(cls, structure_chain, chain_permutation):
        """

        """
        new_coordinates = []
        for chain_id in chain_permutation:
            for atom_coords in structure_chain[chain_id]:
                new_coordinates.append(atom_coords)
        return numpy.array(new_coordinates)

    @classmethod
    def reorderAllCoordinatesByChainLen(cls, structure, selection):
        chains = structure.select(selection).copy().getHierView()

        new_structure = structure.select(selection).copy()
        chain_lens = cls.getChainLengths(new_structure, selection)
        sorted_chain_lens = sorted([(y,x) for x,y in chain_lens.items()])
        coordsets = None
        for len,chid in sorted_chain_lens:
            if coordsets is None:
                coordsets = chains[chid].getCoordsets()
            else:
                coordsets = numpy.concatenate((coordsets, chains[chid].getCoordsets()), axis = 1)
        return coordsets


    @classmethod
    def getChainLengths(cls, structure, selection):
        chains = structure.select(selection).copy().getHierView()
        chain_len_map = {}
        for chain in chains:
            chain_len_map[chain.getChid()] = chain.numAtoms()
        return chain_len_map

    @classmethod
    def calcRMSDMatrix(cls, structure, calculator_type, in_chain_selection ):
        """

        """
        chain_structures = cls.getStructureChains(structure, in_chain_selection) # -> needs padding too!
        chain_ids = chain_structures[0].keys()

        chain_coordsets = structure.select(in_chain_selection).getCoordsets() # -> needs padding

        matrix_data = []

        for i in range(len(chain_structures)-1):
            chain_structure = chain_structures[i]
            min_rmsd = None
            for chain_perm in itertools.permutations(chain_ids):
                new_coords = cls.reorderCoordinates(chain_structure, chain_perm)
                calculator = pyRMSD.RMSDCalculator.RMSDCalculator(calculator_type,
                                                                  numpy.concatenate([[new_coords], chain_coordsets[i+1:]]))
                rmsd = calculator.oneVsFollowing(0)
                if min_rmsd is None:
                    min_rmsd = rmsd
                else:
                    min_rmsd = numpy.minimum(rmsd,min_rmsd)

#                 print "structure:",i, "permutation:",chain_perm, "rmsd:", rmsd
#             print "min rmsd", min_rmsd
            matrix_data.extend(min_rmsd)
        return pyRMSD.condensedMatrix.CondensedMatrix(matrix_data)

