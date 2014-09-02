"""
Created on 27/06/2014

@author: victor
"""
import numpy
import itertools
import pyRMSD.RMSDCalculator
import pyRMSD.condensedMatrix


def combopermutations( elements_list, prefix = []):
    """
    Generator that produces ordered combinations of permutations e.g. for [[1,2],[3,4]] it will produce:
    [[1,2],[3,4]]
    [[2,1],[3,4]]
    [[1,2],[4,3]]
    [[2,1],[4,3]]

    @param elements_list: A list of lists. Each list will be a single permutation.

    @return: a different permutation it time it is called.
    """
    if len(elements_list) <= 1:
        for perm in itertools.permutations(elements_list[0]):
            yield prefix+list(perm)
    else:
        for perm in itertools.permutations(elements_list[0]):
            for perm2 in combopermutations(elements_list[1:], prefix+list(perm)):
                yield  list(perm2)

class ChainMappingRMSDMatrixCalculator:
    """
    Some structures can be formed of repeated units that can be arbitrarily place. This class tries to reorder this units (chains)
    in order to calculate the minimum possible RMSD.

    Ex.1:  Given 3 structures formed of three equal domains:


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

    Ex 2. Different length chains. Suppose we have 4 chains with the following lengths:
    A -> 1
    B -> 2
    C -> 1
    D -> 2

    In this case the permutations to check (for each of the structures) will be:
          1    2
        A, C  B, D
        A, C  D, B
        C, A  B, D
        C, A  D, B
    Note that chains are ordered by length, and that resulting permutations are a combination of the
    possible permutations for groups of chains of the same length (please breathe here).
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
    def getReorderedCoordinatesByLenGroups(cls, structure, selection, len_groups):
        """
        Reordering the chains can affect enough the outcome of the matrix calculation to get different clusterings
        than without ordering. One example is documented in ...
        To diminish changes, all keys are sorted before using them (keys are arbitrary sorted in dics.).
        """
        chains = structure.select(selection).copy().getHierView()
        coordsets = None
        for group_len in sorted(len_groups.keys()):
            for chid in sorted(len_groups[group_len]):
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
    def getChainLenGroups(cls, chain_len_map):
        groups = {}
        for chid in chain_len_map:
            try:
                groups[chain_len_map[chid]].append(chid)
            except:
                groups[chain_len_map[chid]] = [chid]
        return groups

    @classmethod
    def getPermGroups(cls, chain_len_groups):
        perm_groups = []
        for group_len in sorted(chain_len_groups.keys()):
            perm_groups.append(chain_len_groups[group_len])
        return perm_groups

    @classmethod
    def calcRMSDMatrix(cls, structure, calculator_type, in_chain_selection ):
        """

        """
        chain_structures = cls.getStructureChains(structure, in_chain_selection)
        chain_len_map = cls.getChainLengths(structure, in_chain_selection)
        chain_len_groups = cls.getChainLenGroups(chain_len_map)
        chain_coordsets = cls.getReorderedCoordinatesByLenGroups(structure, in_chain_selection, chain_len_groups)
#         perm_file = open("permutations_head_2.txt","w")
        perm_groups = cls.getPermGroups(chain_len_groups)
        matrix_data = []
        for i in range(len(chain_structures)-1):
            chain_structure = chain_structures[i]
            min_rmsd = None
            for chain_perm in combopermutations(perm_groups):
                new_coords = cls.reorderCoordinates(chain_structure, chain_perm)
                calculator = pyRMSD.RMSDCalculator.RMSDCalculator(calculator_type,
                                                                  numpy.concatenate([[new_coords], chain_coordsets[i+1:]]))
                rmsd = calculator.oneVsFollowing(0)
#                 perm_file.write(str(chain_perm)+" "+str(rmsd.tolist())+"\n")
                if min_rmsd is None:
                    min_rmsd = rmsd
                else:
                    min_rmsd = numpy.minimum(rmsd,min_rmsd)
#                 print "structure:",i, "permutation:",chain_perm, "rmsd:", rmsd
#             print "min rmsd", min_rmsd
            matrix_data.extend(min_rmsd)
#         perm_file.close()
        return pyRMSD.condensedMatrix.CondensedMatrix(numpy.array(matrix_data))

