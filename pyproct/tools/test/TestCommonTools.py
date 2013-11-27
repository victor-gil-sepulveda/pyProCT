'''
Created on 16/03/2012

@author: victor
'''
import unittest
import cStringIO
from pyproct.tools.commonTools import merge_files, gen_consecutive_ranges,print_and_flush

class Test(unittest.TestCase):
   
    def test_merge_files(self):
        file_text1 = """Lorem ipsum dolor sit amet, consectetur adipiscing elit.
Proin nec felis mauris, sed semper odio.
Quisque luctus dictum arcu, commodo aliquam lorem porta eget.
"""

        file_text2 = """Aenean facilisis consequat mi, et pretium nibh imperdiet sit amet.
Proin pulvinar eros eu eros molestie a iaculis velit tincidunt.
Curabitur eleifend sem in odio ornare dignissim ut vitae nulla.
"""

        file_text3 = """Suspendisse ac neque magna, nec consequat dui.
Aliquam lacinia lacus a nibh egestas vestibulum.
Pellentesque quis lacus sed lectus porttitor ultricies sit amet nec magna.
Nulla nec justo nunc, ut fringilla sapien.
"""

        file_text4 = """Quisque egestas erat nec nisi viverra hendrerit.
Donec at diam eu ipsum vestibulum convallis ac nec purus.
Curabitur non lorem sed justo commodo placerat.
Donec id arcu odio, in porttitor eros.
Ut ac nibh eu nisi facilisis rhoncus.
Nulla pulvinar mattis lectus, non eleifend libero commodo ut.
"""

        output_text = """Lorem ipsum dolor sit amet, consectetur adipiscing elit.
Proin nec felis mauris, sed semper odio.
Quisque luctus dictum arcu, commodo aliquam lorem porta eget.
Aenean facilisis consequat mi, et pretium nibh imperdiet sit amet.
Proin pulvinar eros eu eros molestie a iaculis velit tincidunt.
Curabitur eleifend sem in odio ornare dignissim ut vitae nulla.
Suspendisse ac neque magna, nec consequat dui.
Aliquam lacinia lacus a nibh egestas vestibulum.
Pellentesque quis lacus sed lectus porttitor ultricies sit amet nec magna.
Nulla nec justo nunc, ut fringilla sapien.
Quisque egestas erat nec nisi viverra hendrerit.
Donec at diam eu ipsum vestibulum convallis ac nec purus.
Curabitur non lorem sed justo commodo placerat.
Donec id arcu odio, in porttitor eros.
Ut ac nibh eu nisi facilisis rhoncus.
Nulla pulvinar mattis lectus, non eleifend libero commodo ut.
"""
        input1 = cStringIO.StringIO(file_text1)
        input2 = cStringIO.StringIO(file_text2)
        input3 = cStringIO.StringIO(file_text3)
        input4 = cStringIO.StringIO(file_text4)
        output = cStringIO.StringIO()
        
        merge_files([input1,input2,input3,input4], output, verbose = False)
        self.assertEqual(output_text,output.getvalue())

    def test_gen_ranges(self):
        expected_range_1 = [0,1,2,3,4]
        expected_range_2 = [5,6,7]
        range_1,range_2 = gen_consecutive_ranges(5,3)
        self.assertItemsEqual(expected_range_1, range_1)
        self.assertItemsEqual(expected_range_2, range_2)
    
    def test_print_and_flush(self):
        handler = cStringIO.StringIO()
        print_and_flush("Hello", handler)
        self.assertEqual(handler.getvalue(), "Hello")
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()