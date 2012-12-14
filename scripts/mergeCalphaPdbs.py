'''
Created on 27/04/2012

@author: victor
'''
import optparse
from pyproclust.tools.pdbTools import create_CA_file
import pyproclust.tools.commonTools as common

if __name__ == '__main__':
    
    parser = optparse.OptionParser(usage='%prog --pdb1 <arg> --pdb2 <arg> -o <arg>',version='1.0')
    parser.add_option('--pdb1', action="store", dest = "pdb1", help="First pdb to be merged",metavar = "my_trajectory_1.pdb")
    parser.add_option('--pdb2', action="store", dest = "pdb2", help="Pdb to be merged at bottom.",metavar = "my_trajectory_2.pdb")
    parser.add_option('-o', action="store", dest = "output", help="Path and name for the output file.",metavar = "out.txt")
    options, args = parser.parse_args()
    
    if not options.pdb1 or not options.pdb2 or not options.output:
        print "Incorrect number of parameters. Please see usage."
    
    ########################################
    # Script Start
    ######################################## 
    
    merged = open(options.output,"w")
    
    ########################################
    # Loading and filtering pdb 1
    ########################################
    common.print_and_flush("Filtering pdb 1...") 
    pdb1 = open(options.pdb1,'r')
    create_CA_file(pdb1,merged)
    pdb1.close()
    common.print_and_flush(" Done\n")
    
    ########################################
    # Loading and filtering pdb 2
    ######################################## 
    common.print_and_flush("Filtering pdb 2...") 
    pdb2 = open(options.pdb2,'r')
    create_CA_file(pdb2,merged)
    pdb2.close()
    common.print_and_flush(" Done\n")
    
    ########################################
    # Finishing script.
    ######################################## 
    merged.close()
