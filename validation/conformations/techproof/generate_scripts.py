'''
Created on 05/06/2014

@author: victor
'''
import pyproct.tools.scriptTools as tools
import os.path

tools.create_directory("scripts/level2")
tools.create_directory("scripts/level1")
tools.create_directory("scripts/level0")

scripts_locations = open("script_locations.txt","w")
all_initial_files = open("files.txt")


BASE_TRAJ_FOLDER = "/gpfs/scratch/bsc72/bsc72476/Victor/2JOF"
BASE_SCRIPT_FOLDER = "scripts"
BASE_CLUSTERING_FOLDER = "compressions"

script_index = 0
all_scripts = []
LEVEL = "level2"
LEVEL2_TEMPLATE = "".join(open(os.path.join(BASE_SCRIPT_FOLDER, "level2_base.json")).readlines())

for file_path in all_initial_files.readlines():
    traj_path = os.path.join(BASE_TRAJ_FOLDER, file_path.strip())
    script_location = os.path.join(BASE_SCRIPT_FOLDER, LEVEL, "script_%s_%d.json"%(LEVEL,script_index))
    workspacelocation = os.path.join(BASE_CLUSTERING_FOLDER,LEVEL,"clustering_%d"%script_index)
    open(script_location,"w").write(LEVEL2_TEMPLATE%(workspacelocation,traj_path))
    all_scripts.append(str(script_location)+"\n")
    script_index += 1

all_initial_files.close()
scripts_locations.close()
open(LEVEL+"_scripts.txt"%all_scripts,"w").write("".join(all_scripts))





