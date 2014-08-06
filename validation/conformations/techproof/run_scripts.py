"""
Created on 05/06/2014

@author: victor
"""

import sys
import os.path
import pyproct.tools.scriptTools as tools

level = sys.argv[1]

execution_range = int(sys.argv[2]) , int(sys.argv[3])+1 #second is included


running_script = """#!/bin/bash
# @ job_name         =  %s
# @ partition        =  bsccv03
# @ initialdir       =  .
# @ output           =  %s.out
# @ error            =  %s.err
# @ total_tasks      =  12
# @ wall_clock_limit =  03:00:00
# @ tasks_per_node   =  12

export PYTHONPATH=/data2/bsc72/clustering/packages:$PYTHONPATH
export LD_LIBRARY_PATH=/data2/apps/OPENMPI/1.6.1/lib:/data2/apps/PYTHON/2.7.5/lib:$LD_LIBRARY_PATH
export PATH=/data2/apps/OPENMPI/1.6.1/bin:$PATH

killall python2.7
mpirun -np 6  /data2/apps/PYTHON/2.7.5/bin/python2.7 -m pyproct.main --mpi %s
"""

exec_path = os.getcwd()

base_folder = {
              "level2":"scripts",
              "level1":"",
              "level0":""
              }

base = os.path.join(base_folder[level], level, "run")
tools.create_directory(base)

for script_number in range(*execution_range):
    script_name = "send_script_%d.run"%script_number
    run_script_file = os.path.join(base, script_name)
    pyproct_script = os.path.join(base_folder[level], level, "script_%s_%d.json"%(level,script_number))
    job_name = "clustering_%s_%d"%(level,script_number)
    open(run_script_file,"w").write(running_script%(job_name,
                                                    os.path.join(base, job_name),
                                                    os.path.join(base, job_name),
                                                    pyproct_script))
    #os.chdir(base)
    #print "send "+run_script_file
    os.system("mnsubmit %s"%run_script_file)
    #os.chdir(exec_path)
