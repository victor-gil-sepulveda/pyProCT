

def removeAllCoordsetsFromStructureLeavingFirst(structure):
    num_conformations = structure.numCoordsets()

    while num_conformations > 1:
        num_conformations-=1
        structure.delCoordset(num_conformations)

def removeAllCoordsetsFromStructure(structure):
    """
    """
    num_conformations = structure.numCoordsets()
    for i in range(num_conformations):
        structure.delCoordset(0)

def getChainPaddedCoordsets(structure, final_len, selection):
    """
    """
    chains = structure.select(selection).copy().getHierView()
    for chain in chains:
        pass


def getMaximumChainLength(structure, selection):
    chains = structure.select(selection).copy().getHierView()
    max_len = 0
    for chain in chains:
        max_len = max(max_len, chain.numAtoms())
    return max_len


