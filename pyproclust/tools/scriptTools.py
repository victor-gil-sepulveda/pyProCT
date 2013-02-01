'''
Created on 19/04/2012

@author: victor
'''
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
    

