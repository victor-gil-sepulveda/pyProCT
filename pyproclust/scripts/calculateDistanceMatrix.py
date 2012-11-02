'''
Created on 27/04/2012

@author: victor
'''
import multiprocessing
import optparse
from pyproclust.matrix.generators.vmdRMSDGenerator import VmdDistanceMatrixGenerator
import pickle
import pyproclust.tools.commonTools as common
import pyproclust.tools.pdbTools as pdb_common

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog --pdb1 <arg> --pdb2 <arg> -o <arg>',version='1.0')
    parser.add_option('--pdb', action="store", dest = "pdb", help="pdb from which we'll calculate the matrix.",metavar = "my_trajectory_1.pdb")
    parser.add_option('--pdb2', action="store", dest = "pdb2", help="pdb from which we'll calculate the matrix.",metavar = "my_trajectory_2.pdb")
    parser.add_option('--text', action="store", dest = "text_matrix_path", help="path for the plain text matrix file. ",metavar = "out.txt")
    parser.add_option('-o', action="store", dest = "output", help="Path and name for the binary matrix file.",metavar = "out.bin")
    parser.add_option('-p', action="store", dest = "number_of_processors", type="int", default = multiprocessing.cpu_count(),help="Number of threads that will be created to process the matrix. By default is the number of cpus.",metavar = "pocessors")
    parser.add_option('--fit_selection', action="store", dest = "fit_selection",  default = "name CA",help="Vmd selection expression for fitting.",metavar = "name CA")
    parser.add_option('--rmsd_selection', action="store", dest = "rmsd_selection", default = "name CA",help="Vmd selection expression for the rmsd calculation.",metavar = "name CA")
    
    options, args = parser.parse_args()
    
    if not options.pdb or not options.output or not options.number_of_processors:
        print "Incorrect number of parameters. Please see usage."
    
    
    ########################################
    # Script Start, matrix calculation
    ######################################## 
    if options.pdb2 == None:
        common.print_and_flush("Creating matrix ...") 
        matrix_generator = VmdDistanceMatrixGenerator(options.pdb, "tmp_matrix", False, options.fit_selection, options.rmsd_selection)
        condensed_matrix = matrix_generator.generate_condensed_matrix(options.number_of_processors,verbose = True)
        common.print_and_flush(" Done\n")
    else:
        common.print_and_flush("Merging trajectories...")
        temporary_merged_trajectory_path =  "merged_tmp.pdb"
        file_handler_out = open(temporary_merged_trajectory_path,"w")
        pdb1_fh = open(options.pdb,'r')
        pdb2_fh = open(options.pdb2,'r')
        pdb_common.create_CA_file (pdb1_fh, file_handler_out)
        pdb_common.create_CA_file (pdb2_fh, file_handler_out)
        file_handler_out.close()
        pdb1_fh.close()
        pdb2_fh.close()
        common.print_and_flush("Done\n")
        
        trajectory_path = temporary_merged_trajectory_path
        common.print_and_flush("Creating matrix ...") 
        matrix_generator = VmdDistanceMatrixGenerator(trajectory_path, "tmp_matrix", False, options.fit_selection, options.rmsd_selection)
        condensed_matrix = matrix_generator.generate_condensed_matrix(options.number_of_processors,verbose = True)
        common.print_and_flush(" Done\n")
    ########################################
    # Writing the binary
    ######################################## 
    output_handler = open(options.output,'w')
    pickle.dump(condensed_matrix,output_handler)
    output_handler.close()
    
    ########################################
    # If required write the plain text matrix
    ######################################## 
    if options.text_matrix_path:
        file_handler = open(options.text_matrix_path,'w')
        condensed_matrix.save(file_handler,low_precision = True)
        file_handler.close()
