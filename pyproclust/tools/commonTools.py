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
    
def convert(my_input):
    """
    Recursively encodes all strings of an input dictionary as UTF-8. Useful to eliminate unicode strings.
    
    @param my_input: A dictionary object.
    
    @return: Encoded dictionary.
    """
    if isinstance(my_input, dict):
        return {convert(key): convert(value) for key, value in my_input.iteritems()}
    elif isinstance(my_input, list):
        return [convert(element) for element in my_input]
    elif isinstance(my_input, unicode):
        return my_input.encode('utf-8')
    else:
        return my_input