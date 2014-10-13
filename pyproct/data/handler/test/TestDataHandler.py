'''
Created on 5/9/2014

@author: victor
'''
import unittest
from pyproct.data.handler.dataHandler import DataHandler

from pyproct.data.handler.test.TestDataLoader import FakeFileLoader

class DataHandlerMock(DataHandler):
    def get_loader(self, data_type):
        return FakeFileLoader

class FakeSourceGenerator():
    def __init__(self, source_list):
        self.source_list = source_list

class TestDataHandler(unittest.TestCase):

    def test_data_handler(self):
        dh = DataHandlerMock({
                          "type":"any", # As our loader is hardcoded it must not check the data type 
                                        # availability
                          "files":[(1,5), (6,9)]
        }, source_generator_class = FakeSourceGenerator)
        
        # We check we have all the elements
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8],
                         list(dh.get_all_elements()))
        
        # Then we check we can get their sources
        self.assertTupleEqual( dh.get_source_of_element(3), (1, 5)) # Element 3 is datum(3) == 4
        self.assertTupleEqual( dh.get_source_of_element(4), (1, 5)) # Element 4 is datum(4) == 5
        self.assertTupleEqual( dh.get_source_of_element(5), (6, 9)) # Element 5 is datum(5) == 6
        self.assertTupleEqual( dh.get_source_of_element(7), (6, 9)) # Element 7 is datum(7) == 8

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()