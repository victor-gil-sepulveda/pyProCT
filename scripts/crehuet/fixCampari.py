'''
Created on 25/03/2014

@author: victor
'''

import os

PELE_MERGER = '/home/victor/workspaces/Python/pyProClust/scripts/pele/peleTrajectoryMerger.py'


def get_files(directory):
    files=os.listdir(directory)
    return sorted(files)

def get_prefix_counts(files, to_this_length):
    prefix = ""
    current_length = len(files)
    i = 1
    prefix_count = {}

    while current_length > to_this_length:
        prefix = files[0][0:i]
        prefix_count[prefix] = 0
        for file_name in files:
            if file_name[0:i] == prefix:
                prefix_count[prefix] += 1
        current_length = prefix_count[prefix]
        i = i + 1
    return prefix_count, prefix

def gen_pdb_folders(to_this_temp):
    pdb_folders = []
    for i in range(to_this_temp+1):
        number = "{:0>2d}".format(i)
        pdb_folders.append("N_0%s"%number)
    return pdb_folders

if __name__ == '__main__':
    for folder in gen_pdb_folders(15):
        prefix_count, last_prefix =  get_prefix_counts(get_files(folder), 1000)

        for i in range(10):
            new_prefix = last_prefix[0:-1]+str(i)
            os.system("python %s -o ../%s.pdb -d %s -p %s --merge_action SEQUENTIAL --append"%(PELE_MERGER,  last_prefix[0:-1], folder, new_prefix))

