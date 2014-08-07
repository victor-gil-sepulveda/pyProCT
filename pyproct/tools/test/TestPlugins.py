'''
Created on 7/8/2014

@author: victor
'''
import unittest
from pyproct.tools.plugins import PluginHandler


test_pckg_exists = False
try:
    import pack1.packmodule
    import pack2.packmodule
    test_pckg_exists = True
except ImportError:
    print "In order to run all tests you need to install the test package 'pack' in pyproct.tools.test.data.entry_point_test"
    
    
class TestPlugins(unittest.TestCase):

    @unittest.skipUnless(test_pckg_exists, "Please install 'pack' in pyproct.tools.test.data.entry_point_test")
    def test_plugin_entry_points(self):
        classes =  PluginHandler.get_classes_from_plugins("testplugin1","pyproct.test.plugin")
        self.assertItemsEqual([x.__name__ for x in classes], ["MatrixTestClass", "AnotherMatrixTestClass"])
        classes =  PluginHandler.get_classes_from_plugins("testplugin2","pyproct.test.plugin")
        self.assertEqual(classes[0].__name__, "AnalysisTestClass")
    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_plugin_entry_points']
    unittest.main()