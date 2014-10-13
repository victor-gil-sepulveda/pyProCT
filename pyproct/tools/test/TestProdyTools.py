"""
Created on 26/06/2014

@author: victor
"""
import unittest
import prody
import numpy
import StringIO
from pyproct.tools.prodyTools import removeAllCoordsetsFromStructure,\
    removeAllCoordsetsFromStructureLeavingFirst
from pyproct.tools.test.data.pdb_data import switched_pdb_data

class Test(unittest.TestCase):
    def test_removeAllCoordsetsFromStructure(self):
        input = StringIO.StringIO(switched_pdb_data)
        pdb_structure = prody.parsePDBStream(input)
        removeAllCoordsetsFromStructure(pdb_structure)
        self.assertEqual(pdb_structure.getCoordsets(),None)

    def test_removeAllCoordsetsFromStructureLeavingFirst(self):
        input = StringIO.StringIO(switched_pdb_data)
        pdb_structure = prody.parsePDBStream(input)
        removeAllCoordsetsFromStructureLeavingFirst(pdb_structure)
        expected = [[[1.0, 2.0, 3.0], [-33.115, 1.294, -1.163]]]
        numpy.testing.assert_array_equal(expected, pdb_structure.getCoordsets())

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()