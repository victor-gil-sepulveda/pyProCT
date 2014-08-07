"""
Created on 25/06/2014

@author: victor
"""
import unittest
import prody
import StringIO
import numpy
from pyproct.driver.handlers.matrix.test.data.pdb_data import pdb_data,\
    chain_padding_proto_3
from pyproct.driver.handlers.matrix.autoChainMappingRMSDMatrixBuilder import ChainMappingRMSDMatrixCalculator,\
    combopermutations

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
        input = StringIO.StringIO(chain_padding_proto_3)
        structure = prody.parsePDBStream(input)
        group_lens = {1: ['A'], 2: ['B', 'D'], 3: ['C', 'E']}
        result =  ChainMappingRMSDMatrixCalculator.getReorderedCoordinatesByLenGroups(structure, "all", group_lens)
        expected =[[[1.0, 2.0, 3.0],
                    [4.0, 5.0, 6.0],
                    [7.0, 8.0, 9.0],
                    [19.0, 20.0, 21.0],
                    [22.0, 23.0, 24.0],
                    [10.0, 11.0, 12.0],
                    [13.0, 14.0, 15.0],
                    [16.0, 17.0, 18.0],
                    [25.0, 26.0, 27.0],
                    [28.0, 29.0, 30.0],
                    [31.0, 32.0, 33.0]],
                   [[41.0, 42.0, 43.0],
                    [44.0, 45.0, 46.0],
                    [47.0, 48.0, 49.0],
                    [59.0, 60.0, 61.0],
                    [62.0, 63.0, 64.0],
                    [50.0, 51.0, 52.0],
                    [53.0, 54.0, 55.0],
                    [56.0, 57.0, 58.0],
                    [65.0, 66.0, 67.0],
                    [68.0, 69.0, 70.0],
                    [71.0, 72.0, 73.0]]]
        numpy.testing.assert_array_equal(expected, result)

    def test_getChainLengths(self):
        input = StringIO.StringIO(chain_padding_proto_3)
        structure = prody.parsePDBStream(input)
        self.assertDictEqual({'A': 1, 'C': 3, 'B': 2, 'E': 3, 'D': 2},
                              ChainMappingRMSDMatrixCalculator.getChainLengths(structure, "all"))

    def test_mult_permutation(self):
        result = [ perm for perm in combopermutations([[1,2], [3,4,5]])]
        expected = [[1, 2, 3, 4, 5], [1, 2, 3, 5, 4], [1, 2, 4, 3, 5], [1, 2, 4, 5, 3], [1, 2, 5, 3, 4], [1, 2, 5, 4, 3], [2, 1, 3, 4, 5], [2, 1, 3, 5, 4], [2, 1, 4, 3, 5], [2, 1, 4, 5, 3], [2, 1, 5, 3, 4], [2, 1, 5, 4, 3]]
        self.assertEqual(result, expected)

    def test_chain_len_groups(self):
        chain_len_map = {'A': 1, 'C': 3, 'B': 2, 'E': 3, 'D': 2}
        expected = {1: ['A'], 2: ['B', 'D'], 3: ['C', 'E']}
        self.assertDictEqual(expected, ChainMappingRMSDMatrixCalculator.getChainLenGroups(chain_len_map))
        self.assertItemsEqual([['A'], ['B', 'D'], ['C', 'E']], ChainMappingRMSDMatrixCalculator.getPermGroups(ChainMappingRMSDMatrixCalculator.getChainLenGroups(chain_len_map)))

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