'''
Created on 16/03/2012

@author: victor
'''

import pyproclust.matrix.generators.generatorTools as common
import pyproclust.tools.pdbTools as pdbtools
import subprocess
from pyproclust.tools.commonTools import merge_files, open_file_handler_list
from pyproclust.matrix.condensedMatrix import load_condensed_matrix

class VmdDistanceMatrixGenerator(object):
    """
    DOCUMENT THIS!!
    """
    def __init__(self,pdb_trajectory_path,matrix_output_file_prefix="partial_out",\
                      symmetry=False, fit_selection="name CA", rmsd_selection="name CA",total_number_of_frames = None):
        """
        TODO: DOCUMENT
        """
        self.pdb_trajectory_path = pdb_trajectory_path
        self.matrix_output_file_prefix = matrix_output_file_prefix
        self.symmetry = symmetry
        self.fit_selection = fit_selection
        self.rmsd_selection = rmsd_selection
        if not total_number_of_frames:
            self.total_frames_of_trajectory = pdbtools.get_number_of_frames(pdb_trajectory_path)
        else:
            self.total_frames_of_trajectory = total_number_of_frames
    
    def generate_condensed_matrix(self, max_num_of_processes, verbose = False):
        """
        Generates a CondensedDistanceMatrix instance!
        """
        tmp_rmsd_file_path = self.matrix_output_file_prefix+"_tmp_rmsd_file"
        files = self.generate_rmsd_files(max_num_of_processes,verbose)
        file_handlers = self.open_file_handler_list(files)
        out_handler = open(tmp_rmsd_file_path,"w")
        
        # Now load matrix
        # It would be great to do this processing in parallel  with a consumer
        # producer scheme :)
        merge_files(file_handlers, out_handler,verbose=verbose)
        out_handler.close()
        out_handler = open(tmp_rmsd_file_path,"r")
        condensed_matrix = load_condensed_matrix(out_handler)
        out_handler.close()
        return condensed_matrix
     
    def gen_vmd_script(self, output, starting_frame, frames_this_iteration):
        """
        Returns a vmd script to calculate a rmsd matrix of the trajectory 'pdb_trajectory'
        with itself. 
        """
        return self.__gen_vmd_script_header(starting_frame,frames_this_iteration,output)\
                        +self.__gen_vmd_script_body()
    
    def generate_rmsd_files(self, processors,verbose = False):
        """
        This function generates a file with the matrix containing the rmsd matrix of a trajectory.
        It creates as much threads as stated in 'processors' and it writes the matrix to the file
        called 'output'. The matrix will be the squared symmmetric matrix if 'symm' is True or
        the upper triangular part if 'symm' is False.
        This is a high level function which is also in charge of distributing the work load but not merging 
        the results. The name of the created files is returned instead, in order to be processed later by the user.
        """
        output_files = []
        processes = []
        ranges = []
        
        ## Get ranges 
        if self.symmetry:
            ranges = common.get_symm_ranges(processors, self.total_frames_of_trajectory)
        else:
            ranges = common.get_nosymm_ranges(processors, self.total_frames_of_trajectory)   
        
        if verbose:
            print "Ranges:",ranges
            
        i = 0
        for r in ranges:
            starting_frame = r[0]
            frames_this_iteration = r[1]
            if verbose:
                print "Launching process from frame", starting_frame
            output_file_name = self.matrix_output_file_prefix+"_"+str(i)
            
            partial_rmsd_matrix_script = self.gen_vmd_script(\
                                            output_file_name,\
                                            starting_frame,\
                                            frames_this_iteration)
            
            processes.append(self.gen_process(i,partial_rmsd_matrix_script))
            
            output_files.append(output_file_name)
            i = i+1
            
        # Wait until all the processes have ended
        if verbose:
            print "All processes were launched. Waiting for termination..."
        for p in processes:
            p.communicate()
            
        # Return the list of files created
        return output_files
    
    def gen_process(self,i,partial_rmsd_matrix_script):
        script_name = self.matrix_output_file_prefix+"_script_"+str(i)+".tcl"
        open(script_name,"w").write(partial_rmsd_matrix_script)
        return subprocess.Popen(["vmd","-dispdev","text","-e",script_name],stdout=subprocess.PIPE)
    
    def open_file_handler_list(self,files):
        return open_file_handler_list(files)

    def __gen_vmd_script_header(self,starting_frame,total_frames_of_this_run,matrix_output_file):
        return """set inpfile "%s"

set outfile  "%s"

set outdescriptor  [open "${outfile}" w]

set starting_frame %d

set total_frames  %d

set sym [lindex $argv 4]

set molid [mol load pdb $inpfile] 

#set nf [molinfo $molid get numframes]

set nf %d

puts "Calculating the RMSD matrix for $nf frames..."
"""%(self.pdb_trajectory_path,matrix_output_file, starting_frame,\
      total_frames_of_this_run, self.total_frames_of_trajectory)
        
    def __gen_vmd_script_body(self):
        """
        Possible selection "name CA"
        """
        script =  """set last_frame [expr $total_frames + $starting_frame]
for { set f1 $starting_frame } { $f1 < $last_frame } { incr f1 1 } {
    
    puts "Calculating frame $f1 of $nf"        
    
    set selection_f1 [atomselect $molid "%s" frame $f1]
    set fit_selection_f1 [atomselect $molid "%s" frame $f1]
    set time_start [clock clicks -milliseconds]
"""%(self.rmsd_selection,self.fit_selection)
        if self.symmetry:
            script = script +"""    set second_start 0
"""
        else:
            script = script +"""    # Set second start to i+1    
    set second_start $f1
    incr second_start 1  
"""
    
        script = script + """    for { set f2 $second_start } { $f2 < $nf } { incr f2 1 } {
        set selection_f2 [atomselect $molid "%s" frame $f2]    
        set fit_selection_f2 [atomselect $molid "%s" frame $f2]
        
        $fit_selection_f2 move [measure fit $selection_f2 $selection_f1]

        puts -nonewline $outdescriptor [format "%%.4f " [measure rmsd $fit_selection_f2 $selection_f1]]

        $selection_f2 delete
        $fit_selection_f2 delete
    }
    puts $outdescriptor ""
    
    puts "It took [expr [clock clicks -milliseconds] - $time_start] ms."

    $selection_f1 delete
    $fit_selection_f1 delete
}
puts "Script finished. $total_frames processed. Closing vmd..."
close $outdescriptor
exit
"""%(self.rmsd_selection,self.fit_selection)
        return script

    

if __name__ == '__main__':
    '''
    import sys
    condensed = VmdDistanceMatrixGenerator(sys.argv[1]).generate_condensed_matrix(1)
    condensed.save(sys.stdout,True)
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
'''