'''
Created on 16/03/2012

@author: victor
'''
import sys

def merge_files(file_handler_list, merged_handler, verbose = True):
    """
    Merges the files which are in the file_handler_list to write them, line per line,
    into the merged_handler in the same order that have in the array.
    
    @param file_handler_list: An array with the file handlers of the files to merge.
    @param merged_handler: The handler of the output file with all the input files merged.
    @param verbose: If true it will print a line every time a file is processed.
    """
    total_files = len(file_handler_list)
    current_file = 1
    if verbose:
        print ""
    for f in file_handler_list:
        if verbose:
            print "Processing file",current_file,"of",total_files
        for line in f:
            merged_handler.write(line)
        current_file = current_file +1 

def gen_consecutive_ranges(num_elems_1,num_elems_2):
    """
    Generates two consecutive ranges from 0 to num_elems_1-1 and from num_elems_1 to
    num_elems_1+num_elems_2 -1 . 
    
    @param num_elems_1: Number of elements of the first range.
    @param num_elems_2: Number of elements of the second range.
    
    @return: The two ranges
    """
    return range(num_elems_1),range(num_elems_1,num_elems_1+num_elems_2)

def print_and_flush(this_string, handler = sys.stdout):
    """
    Prints a string to an opened file handler and makes a flush to ensure it's written.
    
    @param this_string: The string to be written. 
    @param handler: The file handler to write in.
    """
    handler.write(this_string)
    handler.flush()
    
def vararg_callback(option, opt_str, value, parser):
    """
    Parses a list of float numbers. To be used with 'optparse'.
    Got from: http://docs.python.org/2/library/optparse.html
    
    @param option: ...
    @param opt_str: ...
    @param value: ...
    @param parser: ...
    """
    assert value is None
    value = []
    
    def floatable(my_str):
        try:
            float(my_str)
            return True
        except ValueError:
            return False
    
    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(float(arg))
    
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)