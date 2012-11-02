'''
Created on 16/03/2012

@author: victor
'''
import unittest
import pyproclust.tools.pdbTools
import pyproclust.tools.test.data as test_data
import os
import cStringIO


class Test(unittest.TestCase):

    def test_def_get_number_of_frames(self):
        open("test_pdb_for_counting.pdb","w").write(test_data.pdb_1_file_content)
        num_models = pyproclust.tools.pdbTools.get_number_of_frames("test_pdb_for_counting.pdb")
        os.system("rm test_pdb_for_counting.pdb")
        self.assertEqual(num_models,test_data.pdb_1_num_of_models)

    def test_extract_clusters_from_trajectory(self):
        input = cStringIO.StringIO(test_data.pdb_1_file_content)
        output = cStringIO.StringIO()
        frames = [2,5]
        pyproclust.tools.pdbTools.extract_frames_from_trajectory(input, test_data.pdb_1_num_of_models, output, frames,"MODEL","TER")
        self.assertEqual( output.getvalue(),test_data.extracted_pdbs_1)
        
        input = cStringIO.StringIO(test_data.pdb_2_file_content)
        output = cStringIO.StringIO()
        frames = [1,3]
        pyproclust.tools.pdbTools.extract_frames_from_trajectory(input, test_data.pdb_2_num_of_models, output, frames,"MODEL","ENDMDL")
        self.assertEqual( output.getvalue(),test_data.extracted_pdbs_3)
        
        input = cStringIO.StringIO(test_data.pdb_2_file_content)
        output = cStringIO.StringIO()
        frames = [1,3]
        pyproclust.tools.pdbTools.extract_frames_from_trajectory(input, test_data.pdb_2_num_of_models, output, frames,"MODEL","ENDMDL", keep_header = True)
        self.assertEqual( output.getvalue(),test_data.extracted_pdbs_4)
        
    def test_write_a_tfile_model_into_other_tfile(self):
        # We'll write the first and third model
        input_file = cStringIO.StringIO(test_data.pdb_1_file_content)
        output_file = cStringIO.StringIO()
        pyproclust.tools.pdbTools.write_a_tfile_model_into_other_tfile(input_file, output_file, 0,"MODEL","TER", False)
        pyproclust.tools.pdbTools.write_a_tfile_model_into_other_tfile(input_file, output_file, 1,"MODEL","TER", True)
        pyproclust.tools.pdbTools.write_a_tfile_model_into_other_tfile(input_file, output_file, 2,"MODEL","TER", False)
        self.assertEqual( output_file.getvalue(), test_data.extracted_pdbs_2)
    
    def test_read_to_TAG(self):
        input_file = cStringIO.StringIO(test_data.pdb_1_file_content)
        lines = pyproclust.tools.pdbTools.read_to_TAG(input_file,"TER")
        expected_lines = """MODEL 0
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  1.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  1.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  1.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  1.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00            
""".split('\n')
        for i in range(len(lines)):
            self.assertEqual(lines[i][:-1], expected_lines[i])
    
    def test_advance_to_TAG(self):
        input_file = cStringIO.StringIO(test_data.pdb_1_file_content)
        pyproclust.tools.pdbTools.advance_to_TAG(input_file,"MODEL") # MODEL 1
        pyproclust.tools.pdbTools.advance_to_TAG(input_file,"MODEL") # MODEL 2
        pyproclust.tools.pdbTools.advance_to_TAG(input_file,"MODEL") # MODEL 3
        line =  input_file.readline()
        expected_line = "ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  3.00            \n"
        self.assertEqual(line, expected_line)
        
#    def test_advance_to_END_or_TER(self):
#        input_file = cStringIO.StringIO(test_data.pdb_1_file_content)
#        pyproclust.tools.pdbTools.advance_to_END_or_TER(input_file) # to MODEL 1's TER
#        pyproclust.tools.pdbTools.advance_to_END_or_TER(input_file) # to MODEL 2's TER
#        pyproclust.tools.pdbTools.advance_to_END_or_TER(input_file) # to MODEL 3's TER
#        expected_line = "MODEL 4\n"
#        line =  input_file.readline()
#        self.assertEqual(line, expected_line)

    def test_create_CA_file(self):
        input_file = cStringIO.StringIO(test_data.amber_short_contents)
        output_file = cStringIO.StringIO()
        pyproclust.tools.pdbTools.create_CA_file(input_file, output_file)
        self.assertEqual(output_file.getvalue(), test_data.amber_short_ca_contents)
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()