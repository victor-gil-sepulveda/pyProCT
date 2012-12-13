'''
Created on 13/12/2012

@author: victor
'''
import unittest
import pyproclust.tools.plotTools 
import numpy

class Test(unittest.TestCase):

    def test_normalize(self):
        self.fail("TODO")
   
    def test_normalize_in_range(self):
        mylist =[255,40,35,2,245]
        norm = pyproclust.tools.plotTools.normalizeInRange(mylist,0.3,1)
        self.assertItemsEqual(norm,[1.0, 0.40980392156862744, 0.396078431372549, 0.30549019607843136, 0.972549019607843])
        mylist2 = [0.072972369852457405, 0.069252134422401676, 0.076050092873040903, 0.080007843335273554, 0.083781585768440692, 0.071124486539449749, 0.10584964147548799]
        print 1-numpy.array(mylist2)
        norm2 = pyproclust.tools.plotTools.normalizeInRange(1-numpy.array(pyproclust.tools.plotTools.normalizeInRange(mylist2,0,1)),0.3,1)
        print norm2
    
    def test_list2ListWoZeros(self):
        self.fail("TODO")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()