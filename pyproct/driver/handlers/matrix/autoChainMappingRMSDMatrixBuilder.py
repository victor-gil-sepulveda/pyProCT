'''
Created on 27/06/2014

@author: victor
'''
import numpy
import itertools
import pyRMSD.RMSDCalculator
import pyRMSD.condensedMatrix

class ChainMappingRMSDMatrixCalculator:
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
    def calcRMSDMatrix(cls, structure, calculator_type, in_chain_selection ):
        chain_structures = cls.getStructureChains(structure, in_chain_selection)
        chain_ids = chain_structures[0].keys()

        chain_coordsets = structure.select(in_chain_selection).getCoordsets()

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

