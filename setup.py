'''
Created on 25/02/2013

@author: victor
'''
from distutils.core import setup, Extension
import numpy
import distutils.sysconfig

setup(name='pyProClust',
      version='1.0',
      description='',
      author='Victor Alejandro Gil Sepulveda',
      author_email='victor.gil.sepulveda@gmail.com',
      url='https://github.com/victor-gil-sepulveda/pyProClust',
      packages=[
                'pyproct',
                'pyproct.algorithms',
                'pyproct.algorithms.hierarchical',
                'pyproct.algorithms.spectral',
                'pyproct.algorithms.dbscan',
                'pyproct.algorithms.dbscan.cython',
                'pyproct.algorithms.kmedoids',
                'pyproct.algorithms.random',
                'pyproct.algorithms.gromos',
                'pyproct.tools',
                'pyproct.clustering',
                'pyproct.clustering.analysis',
                'pyproct.clustering.comparison',
                'pyproct.clustering.comparison.distrprob',
                'pyproct.clustering.selection',
                'pyproct.clustering.metrics',
                'pyproct.clustering.metrics.cython',
                'pyproct.clustering.filtering',
                'pyproct.protocol',
                'pyproct.protocol.refinement',
                'pyproct.protocol.exploration',
                'pyproct.driver',
                'pyproct.driver.handlers',
                'pyproct.driver.handlers.matrix',
                'pyproct.driver.observer',
                'pyproct.driver.scheduling',
                'pyproct.driver.compressor',
                'pyproct.driver.results',
      ],

      include_dirs = [numpy.get_include(),
                      distutils.sysconfig.get_python_inc()],

      ext_modules=[
                   Extension('pyproct.clustering.metrics.cython.normNCut',[
                                'pyproct/clustering/metrics/cython/normNCut.c'
                   ], extra_compile_args=["-O3","-ffast-math"]),
                   Extension('pyproct.clustering.metrics.cython.boundedCohesion', [
                                'pyproct/clustering/metrics/cython/boundedCohesion.c'
                   ], extra_compile_args=["-O3","-ffast-math"]),
                   Extension('pyproct.clustering.metrics.cython.silhouette',[
                                'pyproct/clustering/metrics/cython/silhouette.c'
                   ], extra_compile_args=["-O3","-ffast-math"]),
                   Extension('pyproct.clustering.metrics.cython.meanMinimumDistance', [
                                'pyproct/clustering/metrics/cython/meanMinimumDistance.c'
                   ], extra_compile_args=["-O3","-ffast-math"]),
                   Extension('pyproct.clustering.metrics.cython.meanMinimumDistance', [
                                'pyproct/clustering/metrics/cython/meanMinimumDistance.c'
                   ], extra_compile_args=["-O3","-ffast-math"]),
                   Extension("pyproct.algorithms.dbscan.cython.cythonDbscanTools", [
                                'pyproct/algorithms/dbscan/cython/cythonDbscanTools.c'
                   ],extra_compile_args=["-O3","-ffast-math"])
      ],

     )
