import sys
import os
import os.path
root_directory = sys.argv[1]

skip_dirs = ["test",".git","build","benchmark"]

def look_for_folders(root,skip_dirs):
    for d in sorted(os.listdir(root)):
        full_path = os.path.join(root,d)
        if os.path.isdir(full_path) and not d in skip_dirs:
            print "'%s',"%full_path.replace(os.sep,".")
            look_for_folders(full_path, skip_dirs)
            
look_for_folders(root_directory,skip_dirs)