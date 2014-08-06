"""
Created on 19/04/2012

@author: victor
"""
import errno
import os.path

def create_directory(directory_path, ensure_writability = False):
    """
    Creates a directory (with subdirs) if it doesn't exist.
    
    @param directory_path: the path of the directory and subdirectories to be created. 
    """
    if ensure_writability:
        if not os.access(os.path.dirname(directory_path), os.W_OK):
            return False
    
    try:
        os.makedirs(directory_path)
        return True
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise
    
    return False

def get_not_repeated_file_name(path_with_file):
    """
    Returns file_name if file_name does not exist. If it exists, it appends an underscore until 
    this new file name does not exist, returning it. For example if "/home/mine/file.txt" exists,
    it will return "/home/mine/_file.txt".
    
    @param path_with_file: complete path with the name of the file.
    """
    directory, file_name = os.path.split(path_with_file)
    file_rename = file_name
    while os.path.exists(directory+"/"+file_rename):
        file_rename = "_"+file_rename
    return directory+"/"+file_rename
    
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
    
