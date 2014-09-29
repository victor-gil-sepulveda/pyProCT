"""
Created on 16/03/2012

@author: victor
"""
import unittest
import pyproct.tools.pdbTools
import pyproct.tools.test.data as test_data
import os
import cStringIO
import numpy
from pyproct.tools.pdbTools import repair_MODEL_ENDMDL_tags,\
    grab_existing_frame_from_trajectory


class Test(unittest.TestCase):

    def test_def_get_number_of_frames(self):
        open("test_pdb_for_counting_1.pdb","w").write(test_data.pdb_1_file_content)
        num_models = pyproct.tools.pdbTools.get_number_of_frames("test_pdb_for_counting_1.pdb")
        os.system("rm test_pdb_for_counting_1.pdb")
        self.assertEqual(num_models,test_data.pdb_1_num_of_models)

    def test_get_remarks_and_filter(self):
        open("test_pdb_for_remarks.pdb","w").write(test_data.pdb_2_file_content)
        
        all_remarks = pyproct.tools.pdbTools.get_remarks("test_pdb_for_remarks.pdb")
        
        self.assertSequenceEqual( all_remarks,
                          [
                               [], 
                               ['REMARK 0 this\n', 'REMARKS 0 is\n', 'REMARKS 0 a remark\n'], 
                               ['REMARKS 0 this\n', 'REMARKS 0 is\n', 'REMARKS 0 a remark\n'], 
                               ['REMARK 400 REGULAR REMARK\n']
                          ])
        
        self.assertSequenceEqual( pyproct.tools.pdbTools.filter_remarks(all_remarks, subset= "ALL"),
                          [
                               [], 
                               ['REMARK 0 this\n', 'REMARKS 0 is\n', 'REMARKS 0 a remark\n'], 
                               ['REMARKS 0 this\n', 'REMARKS 0 is\n', 'REMARKS 0 a remark\n'], 
                               ['REMARK 400 REGULAR REMARK\n']
                          ])
        
        self.assertSequenceEqual( pyproct.tools.pdbTools.filter_remarks(all_remarks, subset= "NOT STANDARD"),
                          [
                               [], 
                               ['REMARK 0 this\n', 'REMARKS 0 is\n', 'REMARKS 0 a remark\n'], 
                               ['REMARKS 0 this\n', 'REMARKS 0 is\n', 'REMARKS 0 a remark\n'], 
                               []
                          ])
        
        self.assertSequenceEqual( pyproct.tools.pdbTools.filter_remarks(all_remarks, subset= "STANDARD"),
                          [
                               [], 
                               [], 
                               [], 
                               ['REMARK 400 REGULAR REMARK\n']
                          ])
        
        self.assertSequenceEqual( pyproct.tools.pdbTools.filter_remarks(all_remarks, subset= "NONE"),
                          [
                               [], 
                               [], 
                               [], 
                               []
                          ])


        os.system("rm test_pdb_for_remarks.pdb")

    def test_get_number_of_atoms(self):
        open("test_pdb_for_counting_2.pdb","w").write(test_data.pdb_1_sub2_file_content)
        num_models = pyproct.tools.pdbTools.get_number_of_atoms("test_pdb_for_counting_2.pdb")
        os.system("rm test_pdb_for_counting_2.pdb")
        self.assertEqual(num_models,test_data.pdb1_num_of_atoms)
        
    def test_extract_frames_from_trajectory_sequentially(self):
        input_file = cStringIO.StringIO(test_data.pdb_1_sub2_file_content)
        output = cStringIO.StringIO()
        frames = [2,5]
        pyproct.tools.pdbTools.extract_frames_from_trajectory_sequentially(input_file, test_data.pdb_1_num_of_models , output, frames,"MODEL","TER")
        self.assertEqual( output.getvalue(),test_data.extracted_pdbs_1)
        
        input_file = cStringIO.StringIO(test_data.pdb_2_file_content)
        output = cStringIO.StringIO()
        frames = [1,3]
        pyproct.tools.pdbTools.extract_frames_from_trajectory_sequentially(input_file, test_data.pdb_2_num_of_models, output, frames,"MODEL","ENDMDL")
        self.assertEqual( output.getvalue(),test_data.extracted_pdbs_3)
        
        input_file = cStringIO.StringIO(test_data.pdb_2_file_content)
        output = cStringIO.StringIO()
        frames = [1,3]
        pyproct.tools.pdbTools.extract_frames_from_trajectory_sequentially(input_file, test_data.pdb_2_num_of_models, output, frames,"MODEL","ENDMDL", keep_header = True)
        self.assertEqual( output.getvalue(),test_data.extracted_pdbs_4)
    
    
    def test_write_a_tfile_model_into_other_tfile(self):
        # We'll write the first and third model
        input_file = cStringIO.StringIO(test_data.pdb_1_sub2_file_content)
        output_file = cStringIO.StringIO()
        pyproct.tools.pdbTools.write_a_tfile_model_into_other_tfile(input_file, output_file, 0,"MODEL","TER", False)
        pyproct.tools.pdbTools.write_a_tfile_model_into_other_tfile(input_file, output_file, 1,"MODEL","TER", True)
        pyproct.tools.pdbTools.write_a_tfile_model_into_other_tfile(input_file, output_file, 2,"MODEL","TER", False)
        self.assertEqual( output_file.getvalue(), test_data.extracted_pdbs_2)
    
    def test_read_to_TAG(self):
        input_file = cStringIO.StringIO(test_data.pdb_1_sub2_file_content)
        lines = pyproct.tools.pdbTools.read_to_TAG(input_file,"TER")
        expected_lines = """MODEL 0
ATOM      3  CA  ILE     3      -0.039   0.638   3.156  1.00  1.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  1.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  1.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  1.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00
""".split('\n')
        for i in range(len(lines)):
            self.assertEqual(lines[i][:-1], expected_lines[i])
    
    def test_advance_to_TAG(self):
        input_file = cStringIO.StringIO(test_data.pdb_1_file_content)
        pyproct.tools.pdbTools.advance_to_TAG(input_file,"MODEL") # MODEL 1
        pyproct.tools.pdbTools.advance_to_TAG(input_file,"MODEL") # MODEL 2
        pyproct.tools.pdbTools.advance_to_TAG(input_file,"MODEL") # MODEL 3
        line =  input_file.readline()
        expected_line = "ATOM      3  CA  ILE     3      -2.039   0.638   3.156  1.00  3.00\n"
        self.assertEqual(line, expected_line)
        
    def test_create_CA_file(self):
        input_file = cStringIO.StringIO(test_data.amber_short_contents)
        output_file = cStringIO.StringIO()
        pyproct.tools.pdbTools.create_CA_file(input_file, output_file)
        self.assertEqual(output_file.getvalue(), test_data.amber_short_ca_contents)
        
    def test_get_model_boundaries(self):
        input_pdb_handler = cStringIO.StringIO(test_data.pdb_1_sub2_file_content)
        boundaries = pyproct.tools.pdbTools.get_model_boundaries(input_pdb_handler)
        expected_boundaries = [[1, 5], [8, 12], [15, 19], [22, 26], [29, 33], [36, 40], [43, 45]]
        numpy.testing.assert_equal(boundaries, expected_boundaries)
        
    def test_repair_MODEL_ENDMDL_tags(self):
        input_pdb_handler = cStringIO.StringIO(test_data.pdb_1_sub2_file_content)
        output_file_handler = cStringIO.StringIO()
        boundaries = [[1, 5], [8, 12], [15, 19], [22, 26], [29, 33], [36, 40], [43, 45]]
        repair_MODEL_ENDMDL_tags(input_pdb_handler, output_file_handler, boundaries)
        self.assertEqual(output_file_handler.getvalue(), test_data.pdb_1_file_content)
        
    def test_grab_existing_frames_from_trajectory(self):
        input_pdb_handler = cStringIO.StringIO(test_data.proto_pdb)
        output_file_handler = cStringIO.StringIO()
        grab_existing_frame_from_trajectory(input_pdb_handler, output_file_handler, 48)
        self.assertEqual(test_data.proto_48_pdb, output_file_handler.getvalue())
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()