'''
Created on 31/01/2012

@author: victor
'''
import unittest
import pyproclust.matrix.generators.test.data as test_data
from pyproclust.matrix.generators.vmdRMSDGenerator import VmdDistanceMatrixGenerator
import cStringIO


class POpenMock(object):
    def __init__(self):
        pass
    def run(self):
        pass
    def communicate(self):
        pass
    
class VmdDistanceMatrixGeneratorMock(VmdDistanceMatrixGenerator):
    def __init__(self,traj_name):
        super(VmdDistanceMatrixGeneratorMock,self).__init__(traj_name,total_number_of_frames=8)
    
    def gen_process(self,i,partial_rmsd_matrix_script):
        return POpenMock()
    
    def open_file_handler_list(self,files):
        handlers = []
        for p in test_data.partial_out_data.all_partials:
            handlers.append(cStringIO.StringIO(p))
        return handlers

class Test(unittest.TestCase):

    def test_gen_vmd_script_w_symmetry(self):
        script_text = test_data.vmd_scripts.symm_script_text
        generator = VmdDistanceMatrixGenerator("clusters.pdb","out",True, "name CA", "name CA",14)
        generated_script = generator.gen_vmd_script("out",0, 14)
        self.assertEqual(script_text,generated_script)

    def test_gen_vmd_script_wout_symmetry(self):
        script_text = test_data.vmd_scripts.no_symm_script_text
        generator = VmdDistanceMatrixGenerator("clusters.pdb","out",False, "name CA", "name CA",14)
        generated_script = generator.gen_vmd_script("out",0, 14)
        self.assertEqual(script_text,generated_script)
        
    def test_generate_condensed_matrix(self):
        #open("tmp_test_data","w").write()
        mock = VmdDistanceMatrixGeneratorMock("lol.pdb")
        condensed_matrix = mock.generate_condensed_matrix(8)
        output = cStringIO.StringIO()
        condensed_matrix.save(output)
        self.assertEqual(test_data.partial_out_data.complete_out_cm_dump,output.getvalue())
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
