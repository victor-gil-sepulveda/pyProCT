"""
Created on 28/05/2013

@author: victor
"""
import unittest
from pyproct.tools.test.data.pdb_data import  merged_1_5, merged_1_5_correlative, proto_pdb
from pyproct.postprocess.actions.representatives import save_representatives


class Test(unittest.TestCase):


    def test_save_representatives(self):
        representatives = [1,5]
        pdbs = ["data/pdb1.pdb", "data/pdb2.pdb", "data/pdb3.pdb"]
        file_path = save_representatives(representatives,
                                         "representatives.pdb",
                                         {"results":"data/","tmp":"data/"},
                                         pdbs,
                                         do_merged_files_have_correlative_models = False,
                                         write_frame_number_instead_of_correlative_model_number = False)
        self.assertEqual("".join(open(file_path,"r").readlines()), merged_1_5_correlative)

        file_path = save_representatives(representatives,
                                         "representatives.pdb",
                                         {"results":"data/","tmp":"data/"},
                                         pdbs,
                                         do_merged_files_have_correlative_models = False,
                                         write_frame_number_instead_of_correlative_model_number = True)
        self.assertEqual("".join(open(file_path,"r").readlines()), merged_1_5)

    def test_extract_prototypes(self):
        self.fail("TODO: Update Test")
        prototypes =  [2,26,48,56,100]
        open("data/prototypes.pdb","w").write(proto_pdb)
        file_path = save_representatives([prototypes[3]],
                                         "prototype.pdb",
                                         {"results":"data/","tmp":"data/"},
                                         ["data/prototypes.pdb"],
                                         do_merged_files_have_correlative_models = False,
                                         write_frame_number_instead_of_correlative_model_number = True)
        print "".join(open(file_path,"r").readlines())

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()