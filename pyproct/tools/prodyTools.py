

def removeAllCoordsetsFromStructureLeavingFirst(structure):
    num_conformations = structure.numCoordsets()

    while num_conformations > 1:
        num_conformations-=1
        structure.delCoordset(num_conformations)

def removeAllCoordsetsFromStructure(structure):
    """
    """
    num_conformations = structure.numCoordsets()
    for _ in range(num_conformations):
        structure.delCoordset(0)



