'''
Created on 25/06/2014

@author: victor
'''

import pyRMSD.RMSDCalculator
import numpy
from pyproct.driver.handlers.trajectoryHandler import TrajectoryHandler
from pyproct.tools.prodyTools import removeAllCoordsetsFromStructure

class ChainMapper:
    """
    Precondition: all models have the same number of chains.
    """
    def __init__(self):
        """

        """
    @classmethod
    def calcChainRMSD(cls, chain_a_coords, chain_b_coords):
        """

        """
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_SERIAL_CALCULATOR", numpy.array([chain_a_coords, chain_b_coords]))

        return calculator.pairwise(0, 1)

    @classmethod
    def twoChainsMap(cls, chain_a, chain_b):
        """
        """
        chain_ids = chain_a.keys()
        number_of_chains = len(chain_ids)
        chain_map = {}
        mapped = []
        for i in range(0,number_of_chains):
            distances = []
            for j in range(0, number_of_chains):
                distances.append((cls.calcChainRMSD(chain_a[chain_ids[i]], chain_b[chain_ids[j]]),j))
            maps_to = min(distances)[1]
            mapped.append(maps_to)
            chain_map[chain_ids[i]] = chain_ids[maps_to]

        valid = (len(set(mapped)) == len(chain_ids))
        if valid:
            return chain_map
        else:
            #return trivial reordering
            trivial_map = {}
            for i in range(0,number_of_chains):
                trivial_map[chain_ids[i]] = chain_ids[i]
            return trivial_map

    @classmethod
    def findMaps(cls, structure_chains):
        """
        """
        maps = []
        for i in range(0,len(structure_chains)):
            maps.append(cls.twoChainsMap(structure_chains[0], structure_chains[i]))
        return maps

    @classmethod
    def getStructureChains(cls, pdb_structure):
        """

        """
        chains = pdb_structure.getHierView()
        number_of_models = pdb_structure.numCoordsets()
        structure_chains = []
        for i in range(number_of_models):
            struct_chain = {}
            for chain in chains:
                struct_chain[chain.getChid()] = chain.getCoordsets(i)
            structure_chains.append(struct_chain)
        return structure_chains

    @classmethod
    def reorderCoordinates(cls, structure_chain, chain_map):
        """
        """
        new_coordinates = []
        for chain_id in sorted(chain_map.keys()):
            for atom_coords in structure_chain[chain_map[chain_id]]:
                new_coordinates.append(atom_coords)
        return numpy.array(new_coordinates)

    @classmethod
    def getRemappedCoordinates(cls, pdb_structure):
        """

        """
        structure_chains = cls.getStructureChains(pdb_structure)
        chain_maps = cls.findMaps(structure_chains)
        new_coordinates = []
        for i in range(len(structure_chains)):
            new_coordinates.append(cls.reorderCoordinates(structure_chains[i], chain_maps[i]))
        return numpy.array(new_coordinates)

class TrajectoryHandlerMapper:

    def __init__(self):
        pass

    def doAutoMapping(self, oldTrajectoryHandler):
        newTrajectoryHandler = TrajectoryHandler(oldTrajectoryHandler.parameters, oldTrajectoryHandler.observer)
        pdb_structure = newTrajectoryHandler.getMergedStructure()
        new_coordinates = ChainMapper.getRemappedCoordinates(pdb_structure)
        removeAllCoordsetsFromStructure(pdb_structure)
        for coordset in new_coordinates:
            pdb_structure.addCoordset(coordset)
        return newTrajectoryHandler
