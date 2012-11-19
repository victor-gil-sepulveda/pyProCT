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

def tile_images_using_montage(file_prefix,num_of_columns,num_of_images,output):
    """
    Uses 'montage' program to generate one image out of some images in a directory
    that start with 'file_prefix'. Images get rearranged in a tiled way.
    """
    num_of_rows =  int(math.ceil(num_of_columns/num_of_images))
    command = "montage -border 0 -geometry 800x -tile %dx%d %s*.png %s"%(num_of_columns,num_of_rows,file_prefix,output)
    os.system(command)

def make_directory(directory_path):
    """
    Creates a directory if it doesn't exist.
    """
    try:
        os.makedirs(directory_path)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise

def save_clusters_as_binary(directory,filename,clustering):
    """
    Saves a clustering as a binary file.
    """
    file_rename = filename
    file_string = directory+"/"+filename
    while  os.path.exists(file_string):
        file_rename = "_"+file_rename
        file_string = directory+"/"+file_rename
    
    file_handler = open(file_string,'w')
    pickle.dump(clustering,file_handler)
    file_handler.close()
    
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
#        print a_file,clusterings[-1].details, len(clusterings[-1].clusters) 
        handler.close()
    return clusterings

#['Spectral', 'DBSCAN', 'GROMOS', 'K-Medoids', 'Random', 'Hierarchical']
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


