"""
Created on 16/03/2012

@author: victor
"""
import sys
import re
import functools

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

def convert_to_utf8(my_input):
    """
    Recursively encodes all strings of an input dictionary as UTF-8. Useful to eliminate unicode strings.

    @param my_input: A dictionary object.

    @return: Encoded dictionary.
    """
    if isinstance(my_input, dict):
        return {convert_to_utf8(key): convert_to_utf8(value) for key, value in my_input.iteritems()}
    elif isinstance(my_input, list):
        return [convert_to_utf8(element) for element in my_input]
    elif isinstance(my_input, unicode):
        return my_input.encode('utf-8')
    else:
        return my_input

def get_parameter_value(key_description, param_dict, default_value):
    """
    Given a parameter description based on dot-separated keys (Ex. data.matrix.type <-> parameters["data"]["matrix"]["type"]) it
    sets the default value for that entry if not defined and returns the default value, or just returns the value if the entry was
    already set. The function fails if any of the intermediate keys does not exist.

    @param key_description:

    @param param_dict:

    @param default_value:

    @return: The value of the dictionary for that key description.
    """
    keys = key_description.split(".")
    tmp_dic = param_dict

    for k in keys[:-1]:
        tmp_dic = tmp_dic[k]

    if not keys[-1] in tmp_dic:
        tmp_dic[keys[-1]] = default_value

    return tmp_dic[keys[-1]]


def remove_comments(string):
    """
    Removes /**/ and // comments from a string (used with the control script).
    From http://stackoverflow.com/questions/2319019/using-regex-to-remove-comments-from-source-files
    """
    string = re.sub(re.compile("/\*.*?\*/",re.DOTALL ) ,"" ,string) # remove all occurance streamed comments (/*COMMENT */) from string
    string = re.sub(re.compile("//.*?\n" ) ,"" ,string) # remove all occurance singleline comments (//COMMENT\n ) from string
    return string

