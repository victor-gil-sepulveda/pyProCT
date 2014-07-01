'''
Created on 25/06/2014

@author: victor
'''
import unittest
import prody
import StringIO
import numpy
from pyproct.driver.handlers.matrix.test.data.pdb_data import pdb_data,\
    chain_padding_proto_3
from pyproct.driver.handlers.matrix.autoChainMappingRMSDMatrixBuilder import ChainMappingRMSDMatrixCalculator
from pyproct.tools.test.data.pdb_data import chain_padding_proto_2

class Test(unittest.TestCase):

    def test_getStructureChains(self):
        expected = [{'A': numpy.array([[ 1.,  2.,  3.]]), 'X': numpy.array([[-33.115,   1.294,  -1.163]])},
                    {'A': numpy.array([[ 4.,  5.,  6.]]), 'X': numpy.array([[-32.555,  -2.5  ,  -5.367]])},
                    {'A': numpy.array([[ 7.,  8.,  9.]]), 'X': numpy.array([[-33.257,   5.28 ,  -8.441]])},
                    {'A': numpy.array([[ 10.,  11.,  12.]]), 'X': numpy.array([[ 32.306,   6.517,  -1.544]])},
                    {'A': numpy.array([[ 13.,  14.,  15.]]), 'X': numpy.array([[ 30.494,  10.39 ,  -3.066]])}]

        input = StringIO.StringIO(pdb_data)
        pdb_structure = prody.parsePDBStream(input)
        result =  ChainMappingRMSDMatrixCalculator.getStructureChains(pdb_structure,"all")
        for i in range(len(expected)):
            numpy.testing.assert_array_equal(expected[i]['A'], result[i]['A'])
            numpy.testing.assert_array_equal(expected[i]['X'], result[i]['X'])

    def test_reorderCoordinates(self):
        chain = {'A': numpy.array([[ 1.,  2.,  3.]]), 'X': numpy.array([[-33.115,   1.294,  -1.163]])}
        chain_map = {"A":"X", "X":"A"}
        numpy.testing.assert_array_equal( [[-33.115,   1.294,  -1.163], [  1., 2., 3.]],
                                          ChainMappingRMSDMatrixCalculator.reorderCoordinates(chain, ["X","A"]))

    def test_reorderAllCoordinatesByChainLen(self):
        input = StringIO.StringIO(chain_padding_proto_2)
        structure = prody.parsePDBStream(input)
        print ChainMappingRMSDMatrixCalculator.reorderAllCoordinatesByChainLen(structure, "all")

    def test_getChainLengths(self):
        input = StringIO.StringIO(chain_padding_proto_3)
        structure = prody.parsePDBStream(input)
        print ChainMappingRMSDMatrixCalculator.getChainLengths(structure, "all")
        print "****\n",ChainMappingRMSDMatrixCalculator.reorderAllCoordinatesByChainLen(structure, "all")

#     def test(self):
#         """
#         Tests if changing the ID changes the coordinates order.
#         """
#         input = StringIO.StringIO(pdb_data)
#         pdb_structure = prody.parsePDBStream(input)
#         chains = pdb_structure.getHierView()
#
#         new_chains = ["X","A"]
#         for i,chain in enumerate(chains):
#             chain.setChid(new_chains[i])
#         output = StringIO.StringIO()
#         prody.writePDBStream(output,pdb_structure)
#         print pdb_structure.getCoordsets()
#         print output.getvalue()

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()