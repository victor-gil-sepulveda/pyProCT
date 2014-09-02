"""
Created on 08/07/2014

@author: victor
"""
import unittest
from pyproct.driver.handlers.matrix.test.data.dihedral_angles_pdb_data import pdb1, expected_dihedrals
import prody
import StringIO
from pyproct.driver.handlers.matrix.dihedralRMSDMatrixCalculator import DihedralRMSDMatrixCalculator,\
    rmsd
import numpy


class Test(unittest.TestCase):

    def test_get_dihedrals(self):
        input = StringIO.StringIO(pdb1)
        pdb_structure = prody.parsePDBStream(input)
        dihedrals =   DihedralRMSDMatrixCalculator.calculateDihedralsForCoordset(pdb_structure,pdb_structure.getCoordsets()[0])
        # We have to get rid off the unknown values!
        numpy.testing.assert_array_almost_equal(numpy.array(expected_dihedrals[1:-1]), numpy.array(dihedrals[2:-2]), 2)

    def test_calc_matrix(self):
        pdb_structure = prody.parsePDB("data/3_models.pdb")
        expected = [ 35.01002624,  47.60315215,  88.64981522,  32.90471145,  87.13023459,  85.76106107]
        product_matrix = DihedralRMSDMatrixCalculator.build(pdb_structure)
#         print "out", product_matrix.get_data()
#         print "out", product_matrix.get_data()
#         print product_matrix.get_data()[0]
#         print product_matrix[0,1]
        numpy.testing.assert_almost_equal(expected, product_matrix.get_data(),8)

    def test_rmsd(self):
        a = [1,2,3,4,5]
        b = [10,20,30,40,50]
        self.assertAlmostEqual(29.8496231132, rmsd(numpy.array(a),numpy.array(b)), 8)
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_build']
    unittest.main()