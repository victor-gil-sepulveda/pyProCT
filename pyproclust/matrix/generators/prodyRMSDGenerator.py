'''
Created on 21/03/2012

@author: victor
'''
from prody.proteins import parsePDB
from prody.ensemble import PDBEnsemble
from multiprocessing import Process,Queue
from pyproclust.matrix.condensedMatrix import CondensedDistanceMatrix
from pyproclust.matrix.generators.generatorTools import get_nosymm_ranges

def calculate_rmsd_in_range(traj_name,shared_data,pdb_data,initial_frame,total_frames_this_run,selection_string):
    trajectory = PDBEnsemble(traj_name)
    trajectory.setAtoms(pdb_data)
    trajectory.addCoordset(pdb_data.getCoordsets()[initial_frame:])
    
    trajectory.select(selection_string) 
    rmsds = []
    i = initial_frame
    for i in range(initial_frame, initial_frame+total_frames_this_run):
        trajectory.setCoords(pdb_data.getCoordsets()[i])
        trajectory.superpose()
        rmsds.extend(trajectory.getRMSDs()[1:])
        trajectory.delCoordset(0)
    
    shared_data.put((i,rmsds))

class ProdyDistanceMatrixGenerator(object):
    '''
    TODOC :D
    '''
    def __init__(self,pdb_trajectory_path, fit_selection = "name CA", rmsd_selection="name CA"):
        """
        Fills everything in.
        """
        self.pdb_trajectory_path = pdb_trajectory_path
        self.fit_selection = fit_selection
        self.rmsd_selection = rmsd_selection
        self.pdb_structure = parsePDB(pdb_trajectory_path)
        self.number_of_frames = len(self.pdb_structure.getCoordsets())
    
    def generate_condensed_matrix(self, max_num_of_processes):
        """
        Returns a condensed matrix as a result of a parallel RMSD calculation with
        'max_num_of_processes' processes.
        """
        ranges = get_nosymm_ranges(max_num_of_processes, self.number_of_frames)
        processes = []
        shared_data = Queue()
         
        # Launch all processes 
        for i in range(len(ranges)):
            traj_name = "trajectory_"+str(i)
            processes.append(self.gen_process(traj_name, shared_data, ranges[i][0], ranges[i][1]))

        # Start everybody
        for p in processes:
            p.start()
        
        # Wait for all to finish
        for p in processes:
            p.join()
        
        # Recreate matrix
        result = []
        while not shared_data.empty():
            result.append(shared_data.get())
        result.sort()
        
        matrix_data = []
        for i in range(len(result)):
            matrix_data.extend(result[i][1])
        condensed_matrix = CondensedDistanceMatrix(matrix_data)
        return condensed_matrix
        

    def gen_process(self, traj_name, shared_data, initial_frame, total_frames_this_run):
        """
        Creates a process to do the hard job of calculating the RMSD (this function has test proxying purpose).
        """
        return Process(target=calculate_rmsd_in_range, args=(traj_name,shared_data,self.pdb_structure, initial_frame,total_frames_this_run,self.rmsd_selection))
    

if __name__ == '__main__':
    import sys
    condensed = ProdyDistanceMatrixGenerator(sys.argv[1]).generate_condensed_matrix(8)
    condensed.save(sys.stdout,True)
    """
    for amber_short.pdb    
    0.6433 0.8617 0.9032 1.0380 0.8893 0.6157 0.6869 0.8225 1.0331 1.0892 1.0714 
    0.8412 0.8098 1.0032 0.8385 0.5830 0.7427 0.9783 0.8659 1.0352 1.0417 
    1.0999 1.4431 0.9664 0.8423 1.0616 1.1735 1.2850 1.0537 1.2200 
    0.8226 1.0650 0.7687 0.8181 0.8753 0.7983 0.9643 0.8349 
    1.2578 1.0073 0.9528 1.0867 0.9220 1.3512 1.1609 
    0.7453 1.1024 1.2526 1.0756 1.0810 1.0502 
    0.7823 0.9234 0.8889 0.9568 0.9760 
    0.7883 1.0015 1.1214 1.0354 
    1.0227 0.9876 1.0371 
    0.9208 0.8179 
    0.7175
    """