'''
Created on 13/03/2012

@author: victor
'''
import unittest
from pyproclust.matrix.matrixCommon import generate_row


class Test(unittest.TestCase):

    def test_generate_row(self):
        row = []
        generate_row(" 1 2 3 4 5  \t 6", row)
        self.assertEqual( row, [1,2,3,4,5,6] )

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()