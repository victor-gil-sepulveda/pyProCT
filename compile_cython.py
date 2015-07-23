"""
Created on 23/7/2015

@author: victor
"""

import os

installers = [
"pyproct/clustering/algorithms/spectral/cython/",
"pyproct/clustering/algorithms/dbscan/cython/",
"pyproct/clustering/evaluation/metrics/cython/"]

current_dir = os.getcwd()

for inst_path in installers:
    os.chdir(inst_path)
    os.system("python setup.py build_ext --inplace")

os.chdir(current_dir)