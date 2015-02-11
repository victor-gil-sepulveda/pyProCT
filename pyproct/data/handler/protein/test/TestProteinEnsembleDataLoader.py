'''
Created on 5/9/2014

@author: victor
'''

import os
import unittest
from pyproct.data.handler.dataSource import DataSource
from pyproct.driver.parameters import ProtocolParameters
import pyproct.data.handler.protein.test.data as test_data
from pyproct.data.handler.protein.proteinEnsembleDataLoader import ProteinEnsembleDataLoader

class TestProteinEnsembleDataLoader(unittest.TestCase):
    
    def setup_loader(self):
        """
        Builds a loader with 2 preloaded pdbs from the data folder.
        """
        loader = ProteinEnsembleDataLoader(ProtocolParameters({"matrix":{}}))
        source1 = DataSource(os.path.join(test_data.__path__[0], "pdb1.pdb"))
        loader.load(source1)
        source2 = DataSource({
                              "source":os.path.join(test_data.__path__[0], "pdb2.dcd"),
                              "atoms_source": os.path.join(test_data.__path__[0], "pdb2.pdb")
        })
        loader.load(source2)
        return loader

    def test_load(self):
        """
        If we already loaded 2 pdbs of 2 conformations, the elements corresponfding to 
        the next loaded pdb of 2 confs, 2 atoms will be [4,5] (elements start in 0).
        It indirectly tests the dcd loading (which would be a responsabilit of the
        data test indeed).
        """
        loader = self.setup_loader()
        source3 = DataSource(os.path.join(test_data.__path__[0], "pdb3.pdb"))
        final_range = loader.load(source3)
        self.assertEqual([4,5], list(final_range))
        self.assertEqual(source3.source['number_of_conformations'], 2)
        self.assertEqual(source3.source['number_of_atoms'], 2)
        
    def test_load_with_base_selection(self):
        """
        Fine grain case of test_load using base_selection in the file descriptor.
        """
        loader = self.setup_loader()
        source3 = DataSource({
                              "source": os.path.join(test_data.__path__[0], "pdb3.pdb"),
                              "base_selection": "resnum 3"
        })
        final_range = loader.load(source3)
        self.assertSequenceEqual([4,5], list(final_range))
        self.assertEqual(source3.source['number_of_conformations'], 2)
        self.assertEqual(source3.source['number_of_atoms'], 1)
        self.assertEqual(loader.loaded_data[-1].getResnames(), ['ILE'])
        
        # Closing this loader must raise a prody exception (different number of atoms
        # for some conformations)
        self.assertRaises(ValueError, loader.close)

    def get_x_coords(self, coordsets):
        x_coords = []
        for conf in coordsets:
            for atom in conf:
                x_coords.append( int(atom[0]))
        return x_coords
    
    def test_close(self):
        """
        Tests that the merged trajectory has been created.
        The coordinates of pdbs 1 to 3 all start with an integer from 
        1 to 8, which is what we are going to test. 
        """
        loader = self.setup_loader()
        data = loader.close()
        coordinates = data.getCoordinates()
        x_coords = self.get_x_coords(coordinates)
                
        self.assertSequenceEqual( range(1,9), x_coords)
    
    def test_close_with_model_remark(self):
        """
        Tests that model number and remark info was acquired
        """
        loader = self.setup_loader()
        source3 = DataSource(os.path.join(test_data.__path__[0], "pdb3.pdb"))
        loader.load(source3)
        loader.close()
        expected_models = [1,2,1,2,8,9] #[1,2,5,6,8,9] -> this if whe have a pdb
        expected_remarks = [
                            ["REMARK 400 HELLO\n"], 
                            ["REMARK 400 WORLD\n"],
                            [],
                            [],
                            [],
                            ["REMARKS NINETH MODEL\n"]
                           ]
        self.assertSequenceEqual(expected_remarks, loader.model_remarks)
        self.assertSequenceEqual(expected_models, loader.model_numbers)
    
    def test_data(self):
        """
        Tests that the mandatory method 'get_element' works as expected.
        """
        loader = self.setup_loader()
        data = loader.close()
        x_coords = self.get_x_coords(data.get_element(2).getCoordsets())
        self.assertSequenceEqual( [5,6], x_coords)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()