'''
Created on 5/9/2014

@author: victor
'''
import unittest
from pyproct.data.handler.elementRange import ElementRange
from pyproct.data.handler.dataLoader import DataLoader


class FakeFileLoader(DataLoader):
    LOADER_TYPE = "test::stuff"
    
    def load_data_from_source(self, data_source):
        e_range = ElementRange(data_source[0], data_source[1])
        return e_range, len(e_range)
    
    def close(self):
        super(FakeFileLoader,self).close()
        inflated_list = []
        [inflated_list.extend(list(e_range)) for e_range in self.loaded_data]
        return inflated_list


class Test(unittest.TestCase):

    def test_all(self):
        ffl = FakeFileLoader(None)
        # Init worked
        self.assertTrue(hasattr(ffl, "loaded_data"))

        # If the loader is empty, closing it triggers exit()
        with self.assertRaises(SystemExit):
            ffl.close()
        
        # If we load 2 datum, we have to get the correct ranges
        ffl.load((1,5))
        e_range = ffl.load((6,9))
        self.assertEqual([5, 6, 7, 8], list(e_range))
        
        # If we close after getting some elements, we have the sum
        # of all of them (this tests FakeFileLoader to complement 
        # the test)
        merged_data = ffl.close()
        self.assertEqual([1, 2, 3, 4, 5, 6, 7, 8, 9], merged_data)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()