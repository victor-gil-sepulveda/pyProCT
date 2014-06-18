'''
Created on 05/06/2014

@author: victor
'''

import sys
import os.path
import pyproct.tools.scriptTools as tools

level = sys.argv[1]

execution_range = int(sys.argv[2]) , int(sys.argv[3])+1 #second is included

exec_path = os.getcwd()

base_folder = {
              "level2":"scripts",
              "level1":"scripts",
              "level0":""
              }

base = os.path.join(base_folder[level], level, "run")
tools.create_directory(base)

for script_number in range(*execution_range):
    pyproct_script = os.path.join(base_folder[level], level, "script_%s_%d.json"%(level,script_number))
    os.system("python -m pyproct.main %s"%pyproct_script)
