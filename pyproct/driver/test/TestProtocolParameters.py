"""
Created on 15/07/2014

@author: victor
"""
import unittest
from pyproct.driver.parameters import ProtocolParameters


class Test(unittest.TestCase):


    def test_to_dic(self):
        input = ProtocolParameters({
                                   "ciao":ProtocolParameters({
                                                              "come": "ti",
                                                              "chiami": 4}),
                                   "io_mi":{
                                            "chiamo":87,
                                            "Victor":[1,2,3,4]
                                            }
                                   })
        expected = {
                    'ciao': {
                             'come': 'ti',
                             'chiami': 4
                             },
                    'io_mi': {
                              'chiamo': 87,
                              'Victor': [1, 2, 3, 4]
                              }
                    }
        self.assertDictEqual(expected, ProtocolParameters.to_dict(input))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_to_dic']
    unittest.main()