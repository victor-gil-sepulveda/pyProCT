'''
Created on 23/9/2014

@author: victor
'''
import unittest
import prody
from pyproct.data.handler.protein.proteinEnsembleData import ProteinEnsembleData
import pyproct.data.handler.protein.test.data as test_data
import os

class Test(unittest.TestCase):

    def test_contents(self):
        
        ped = ProteinEnsembleData(prody.parsePDB(os.path.join(test_data.__path__[0],"pdb123.pdb")), {})
        self.assertEqual(ped.get_num_elements(), 6)
        
        def get_x_coords(coordsets):
            x_coords = []
            for conf in coordsets:
                for atom in conf:
                    x_coords.append(int(atom[0]))
            return x_coords
        
        self.assertItemsEqual(range(1,13), get_x_coords(ped.getCoordinates()))
        
        self.assertItemsEqual([3, 4, 9, 10, 11, 12], 
                              get_x_coords(ped.get_elements([1,4,5]).getCoordsets()))
        
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()