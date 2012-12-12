'''
Created on 19/04/2012

@author: victor
'''
import pickle
import os.path
import errno
import math

"""
Some functions without a good place to be put :)
"""

def get_directories_list(base_path):
    directories = []
    for dirname, dirnames, filenames in os.walk(base_path,topdown=False):
        try:
            directories.append(dirname.split("/")[1])
        except:
            pass    
    return list(set(directories))


def make_directory(directory_path):
    """
    Creates a directory if it doesn't exist.
    """
    try:
        os.makedirs(directory_path)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise

# def save_clusters_as_binary(directory,filename,clustering):
#     """
#     Saves a clustering as a binary file.
#     """
#     file_rename = filename
#     file_string = directory+"/"+filename
#     while  os.path.exists(file_string):
#         file_rename = "_"+file_rename
#         file_string = directory+"/"+file_rename
#     
#     file_handler = open(file_string,'w')
#     pickle.dump(clustering,file_handler)
#     file_handler.close()

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
    
def load_binary_clusters(clusterings_dir):
    clusterings = []
    clustering_files = []
    
    files = os.listdir(clusterings_dir) 
    for filename in files:
        if ".bin" in filename:
            clustering_files.append(clusterings_dir+"/"+filename)
    clustering_files.sort()
    for a_file in clustering_files:
        handler = open(a_file, "r")
        clusterings.append(pickle.load(handler))
        handler.close()
    return clusterings

def classify_generated_clusters(tags,clusterings):
    counter = {}
    for t in tags:
        counter[t] = 0
    for clustering in clusterings:
        for t in tags:
            if t in clustering.details:
                counter[t] += 1
    return counter, len(clusterings)

def save_report(filename,analyzer):
    """
    Saves the clustering report to a file.
    """
    make_directory('./reports')
    file_handler = open('./reports'+"/"+filename,'w')
    file_handler.write(analyzer.generate_report())
    file_handler.close()

# def tile_images_using_montage(file_prefix,num_of_columns,num_of_images,output):
#     """
#     Uses 'montage' program to generate one image out of some images in a directory
#     that start with 'file_prefix'. Images get rearranged in a tiled way.
#     """
#     num_of_rows =  int(math.ceil(num_of_columns/num_of_images))
#     command = "montage -border 0 -geometry 800x -tile %dx%d %s*.png %s"%(num_of_columns,num_of_rows,file_prefix,output)
#     os.system(command)

