'''
Created on 25/06/2014

@author: victor
'''
import unittest
from pyproct.preprocessing.test.data.pdb_data import pdb_data, switched_pdb_data
import prody
import StringIO
from pyproct.preprocessing.chainMapping import ChainMapper
import numpy
from math import sqrt
import types

class Test(unittest.TestCase):

    def test_getStructureChains(self):
        expected = [{'A': numpy.array([[ 1.,  2.,  3.]]), 'X': numpy.array([[-33.115,   1.294,  -1.163]])},
                    {'A': numpy.array([[ 4.,  5.,  6.]]), 'X': numpy.array([[-32.555,  -2.5  ,  -5.367]])},
                    {'A': numpy.array([[ 7.,  8.,  9.]]), 'X': numpy.array([[-33.257,   5.28 ,  -8.441]])},
                    {'A': numpy.array([[ 10.,  11.,  12.]]), 'X': numpy.array([[ 32.306,   6.517,  -1.544]])},
                    {'A': numpy.array([[ 13.,  14.,  15.]]), 'X': numpy.array([[ 30.494,  10.39 ,  -3.066]])}]

        input = StringIO.StringIO(pdb_data)
        pdb_structure = prody.parsePDBStream(input)
        result =  ChainMapper.getStructureChains(pdb_structure)
        for i in range(len(expected)):
            numpy.testing.assert_array_equal(expected[i]['A'], result[i]['A'])
            numpy.testing.assert_array_equal(expected[i]['X'], result[i]['X'])

    def test_twoChainsMap(self):
        chainMapperStub = ChainMapper()
        def new_calcChainRMSD(self, chain_a, chain_b):
            return sqrt(((chain_a[0] - chain_b[0])**2).sum()) # euclidean distance
        setattr(ChainMapper, "calcChainRMSD", types.MethodType(new_calcChainRMSD, ChainMapper))

        chains =[   {'A': numpy.array([[ 1.,  2.,  3.]]), 'X': numpy.array([[-33.115,   1.294,  -1.163]])},
                    {'A': numpy.array([[ 4.,  5.,  6.]]), 'X': numpy.array([[-32.555,  -2.5  ,  -5.367]])},
                    {'A': numpy.array([[-33.257,   5.28 ,  -8.441]]), 'X': numpy.array([[ 7.,  8.,  9.]])}, # Switched chains
                    {'A': numpy.array([[ 7.,  8.,  9.]]), 'X': numpy.array([[ 7.,  8.,  9.]])},# Must return None
                    ]

        self.assertDictEqual( {'A': 'A', 'X': 'X'}, chainMapperStub.twoChainsMap(chains[0], chains[1]))
        self.assertDictEqual( {'A': 'X', 'X': 'A'}, chainMapperStub.twoChainsMap(chains[0], chains[2]))
        self.assertDictEqual( {'A': 'A', 'X': 'X'}, chainMapperStub.twoChainsMap(chains[0], chains[3]))

    def test_findMaps(self):
        chainMapperStub = ChainMapper()
        def new_calcChainRMSD(self, chain_a, chain_b):
            return sqrt(((chain_a[0] - chain_b[0])**2).sum()) # euclidean distance
        setattr(ChainMapper, "calcChainRMSD", types.MethodType(new_calcChainRMSD, ChainMapper))
        chains =[   {'A': numpy.array([[ 1.,  2.,  3.]]), 'X': numpy.array([[-33.115,   1.294,  -1.163]])},
                    {'A': numpy.array([[ 4.,  5.,  6.]]), 'X': numpy.array([[-32.555,  -2.5  ,  -5.367]])},
                    {'A': numpy.array([[-33.257,   5.28 ,  -8.441]]), 'X': numpy.array([[ 7.,  8.,  9.]])}, # Switched chains
                    {'A': numpy.array([[ 7.,  8.,  9.]]), 'X': numpy.array([[ 7.,  8.,  9.]])},# Must return None
                    ]
        expected = [{'A': 'A', 'X': 'X'}, {'A': 'A', 'X': 'X'}, {'A': 'X', 'X': 'A'}, {'A': 'A', 'X': 'X'}]
        for result in  chainMapperStub.findMaps(chains):
            self.assertIn(result, expected)
        self.assertEqual(4, len(chainMapperStub.findMaps(chains)))

    def test_reorderCoordinates(self):
        chain = {'A': numpy.array([[ 1.,  2.,  3.]]), 'X': numpy.array([[-33.115,   1.294,  -1.163]])}
        chain_map = {"A":"X", "X":"A"}
        numpy.testing.assert_array_equal( [[-33.115,   1.294,  -1.163], [  1., 2., 3.]], ChainMapper.reorderCoordinates(chain, chain_map))

    def test_getRmappedCoordinates(self):
        input = StringIO.StringIO(switched_pdb_data)
        switched_pdb_structure = prody.parsePDBStream(input)

        chainMapperStub = ChainMapper()
        def new_calcChainRMSD(self, chain_a, chain_b):
            return sqrt(((chain_a[0] - chain_b[0])**2).sum()) # euclidean distance
        setattr(ChainMapper, "calcChainRMSD", types.MethodType(new_calcChainRMSD, ChainMapper))
        expected = [[[1.0, 2.0, 3.0], [-33.115, 1.294, -1.163]],
                    [[4.0, 5.0, 6.0], [-32.555, -2.5, -5.367]],
                    [[7.0, 8.0, 9.0], [-33.257, 5.28, -8.441]],
                    [[10.0, 11.0, 12.0], [-32.306, 6.517, -1.544]],
                    [[13.0, 14.0, 15.0], [-30.494, 10.39, -3.066]]]

        numpy.testing.assert_array_equal(expected, chainMapperStub.getRemappedCoordinates(switched_pdb_structure))

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