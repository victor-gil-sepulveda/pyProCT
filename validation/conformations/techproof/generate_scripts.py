"""
Created on 05/06/2014

@author: victor
"""
import pyproct.tools.scriptTools as tools
import os.path

import sys

LEVEL = sys.argv[1]

tools.create_directory("scripts/level2")
tools.create_directory("scripts/level1")
tools.create_directory("scripts/level0")

all_initial_files = open(sys.argv[2])
print "Working with filelist: %s"%sys.argv[2]

BASE_TRAJ_FOLDER = "/gpfs/scratch/bsc72/bsc72476/Victor/2JOF"
BASE_SCRIPT_FOLDER = "scripts"
BASE_CLUSTERING_FOLDER = "compressions"

script_index = 0
all_scripts = []


LEVEL_TEMPLATE = "".join(open(os.path.join(BASE_SCRIPT_FOLDER, "%s_base.json"%LEVEL)).readlines())

for file_path in all_initial_files.readlines():
    traj_path = os.path.join(BASE_TRAJ_FOLDER, file_path.strip())
    script_location = os.path.join(BASE_SCRIPT_FOLDER, LEVEL, "script_%s_%d.json"%(LEVEL,script_index))
    workspacelocation = os.path.join(BASE_CLUSTERING_FOLDER,LEVEL,"clustering_%d"%script_index)
    open(script_location,"w").write(LEVEL_TEMPLATE%(workspacelocation,traj_path))
    all_scripts.append(str(script_location)+"\n")
    script_index += 1

all_initial_files.close()
open(LEVEL+"_scripts.txt"%all_scripts,"w").write("".join(all_scripts))





