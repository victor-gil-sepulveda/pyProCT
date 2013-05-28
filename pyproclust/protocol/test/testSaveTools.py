'''
Created on 28/05/2013

@author: victor
'''
import unittest
from pyproclust.tools.test.data.pdb_data import premerged_pdb_1, premerged_pdb_2,\
    premerged_pdb_3, merged_pdb, merged_1_5, merged_renumbered_pdb,\
    merged_1_5_correlative, proto_pdb
from pyproclust.protocol.saveTools import merge_pdbs, save_representatives
from pyproclust.driver.handlers.trajectoryHandler import TrajectoryHandler


class Test(unittest.TestCase):

    def test_merge_pdbs(self):
        pdbs = {"data/pdb1.pdb":premerged_pdb_1,
                "data/pdb2.pdb":premerged_pdb_2,
                "data/pdb3.pdb":premerged_pdb_3} 
        
        for pdb_file in pdbs:
            open(pdb_file,"w").write(pdbs[pdb_file])
        
        merge_pdbs(sorted(pdbs.keys()), 
                   "data/merged.pdb", 
                   do_merged_files_have_correlative_models=False)
        merged = "".join(open("data/merged.pdb","r").readlines())
        self.assertEqual(merged, merged_pdb)
        
        trajectoryHandler = TrajectoryHandler({"pdbs":["data/pdb1.pdb", "data/pdb2.pdb", "data/pdb3.pdb"]}, 
                                              None)
        merge_pdbs(trajectoryHandler,
                   "data/merged.pdb",
                   do_merged_files_have_correlative_models=True)
        merged = "".join(open("data/merged.pdb","r").readlines())
        self.assertEqual(merged, merged_renumbered_pdb)
        
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