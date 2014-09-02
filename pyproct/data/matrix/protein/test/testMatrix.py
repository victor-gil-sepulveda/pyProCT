"""
Created on 05/02/2014

@author: victor
"""
import unittest
from pyproct.driver.handlers.matrix.rmsdMatrixBuilder import RMSDMatrixBuilder
import numpy

class TrajectoryHandlerStub():
    def __init__(self, sel_result_pairs ):
        self.results = {}
        for key, result in sel_result_pairs:
            self.results[key] = result

    def getSelection(self, selection):
        return self.results[selection]

class Test(unittest.TestCase):

    def test_select_one_atom(self):
        traj_handler = TrajectoryHandlerStub([("one atom",[[[1,2,3]]]),
                                             ("another atom", [[[4,5,6]]]),
                                             ("two atoms", [[[4,5,6],[7,8,9]]]),
                                             ("no atoms", [[]])])

        self.assertItemsEqual([1,2,3], RMSDMatrixBuilder.select_one_atom(traj_handler, "one atom"))

        with self.assertRaises(SystemExit):
            RMSDMatrixBuilder.select_one_atom(traj_handler, "two atoms")

        with self.assertRaises(SystemExit):
            RMSDMatrixBuilder.select_one_atom(traj_handler, "no atoms")

    def test_process_group(self):
        matrix_parameters = {
                             "symmetries":{
                                           "first_equivalence":{
                                                                "common": "my common text",
                                                                "equivalences": [
                                                                                 ["one atom","another atom"],
                                                                                 ["yet another atom","final atom"]
                                                                                 ]
                                                                },
                                           "second_equivalence":{
                                                                "common": "my common text",
                                                                "equivalences": [
                                                                                 ["one atom","yet another atom"]
                                                                                 ]
                                                                }
                              }
                             }

        traj_handler = TrajectoryHandlerStub([   ("my common text and one atom",numpy.array([[[1,2,3]]])),
                                                 ("my common text and another atom", numpy.array([[[4,5,6]]])),
                                                 ("my common text and yet another atom", numpy.array([[[7,8,9]]])),
                                                 ("my common text and final atom", numpy.array([[[10,11,12]]]))])

        calc_selection_coordsets = numpy.array([
                                                [
                                                 [0,0,0],
                                                 [1,2,3],
                                                 [-1,-2,-3],
                                                 [4,5,6],
                                                 [-4,-5,-6],
                                                 [7,8,9],
                                                 [-7,-8,-9],
                                                 [10,11,12],
                                                 [-10,-11,-12]
                                              ]
                                            ])


        self.assertItemsEqual([(1, 3), (5, 7)],RMSDMatrixBuilder.process_group("first_equivalence",
                                              matrix_parameters,
                                              traj_handler,
                                              calc_selection_coordsets))

        self.assertItemsEqual([[(1, 3), (5, 7)],[(1,5)]] ,RMSDMatrixBuilder.process_symm_groups(matrix_parameters,
                                              traj_handler,
                                              calc_selection_coordsets))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()