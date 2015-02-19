"""
Created on 26/07/2012

@author: victor
"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

if __name__ == '__main__': # Comp. with sphynx
    setup(
        cmdclass = {'build_ext': build_ext},
        include_dirs = [np.get_include()],
        ext_modules = [
                       Extension("cohesion", ["cohesion.pyx"], extra_compile_args=["-O3","-ffast-math"]),
                       Extension("silhouette", ["silhouette.pyx"], extra_compile_args=["-O3","-ffast-math"]),
                       Extension("graph.tools", ["graph/tools.pyx"], extra_compile_args=["-O3","-ffast-math"]),
                       Extension("graph.ratioCut", ["graph/ratioCut.pyx"], extra_compile_args=["-O3","-ffast-math"]),
                       Extension("graph.minMaxCut", ["graph/minMaxCut.pyx"], extra_compile_args=["-O3","-ffast-math"]),
                       Extension("graph.nCut", ["graph/nCut.pyx"], extra_compile_args=["-O3","-ffast-math"])
        ])
