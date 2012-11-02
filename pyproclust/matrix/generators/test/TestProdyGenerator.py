'''
Created on 21/03/2012

@author: victor
'''
import unittest
from pyproclust.matrix.generators.prodyRMSDGenerator import ProdyDistanceMatrixGenerator
from multiprocessing.process import Process
import cStringIO



def parsePDB(filename):
    return ProdyPDBStructureMock(filename)

class ProdyPDBStructureMock(object):
    
    def __init__ (self,filename):
        self.filename = filename
    
    def getCoordsets(self):
        return [[0]*3*10]*10

def do_something_with_matrix(traj_name,shared_data,pdb_data,initial_frame,total_frames_this_run,selection_string):
    total_frames = len(pdb_data.getCoordsets())
    rmsds = []
    for i in range(initial_frame, initial_frame+total_frames_this_run):
        rmsds.extend([i]*(total_frames-i-1))
    shared_data.put((initial_frame,rmsds))

class ProdyDistanceMatrixGeneratorMock(ProdyDistanceMatrixGenerator):
    
    def __init__(self,pdb_trajectory_path):
        self.pdb_trajectory_path = pdb_trajectory_path
        self.fit_selection = "name CA"
        self.rmsd_selection = "name CA"
        self.pdb_structure = parsePDB(pdb_trajectory_path)
        self.number_of_frames = len(self.pdb_structure.getCoordsets())
        
    def gen_process(self, traj_name, shared_data, initial_frame, total_frames_this_run):
        return Process(target=do_something_with_matrix, args=(traj_name,shared_data,self.pdb_structure, initial_frame,total_frames_this_run,self.rmsd_selection))
    
class Test(unittest.TestCase):

    def test_generate_condensed_matrix(self):
        expected_matrix = """0 0 0 0 0 0 0 0 0 
1 1 1 1 1 1 1 1 
2 2 2 2 2 2 2 
3 3 3 3 3 3 
4 4 4 4 4 
5 5 5 5 
6 6 6 
7 7 
8 
"""
        output = cStringIO.StringIO()
        p = ProdyDistanceMatrixGeneratorMock("trajectory.pdb")
        condensed = p.generate_condensed_matrix(8)
        condensed.save(output)
        self.assertEquals(expected_matrix,output.getvalue())
        
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_generate_condensed_matrix']
    unittest.main()