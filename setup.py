'''
Created on 25/02/2013

@author: victor
'''
from distutils.core import setup, Extension

setup(name='pyProClust',
      version='1.0',
      description='',
      author='Victor Alejandro Gil Sepulveda',
      author_email='victor.gil.sepulveda@gmail.com',
      url='https://github.com/victor-gil-sepulveda/pyProClust',
      packages=[
                'pyproclust',
                'pyproclust.algorithms',
                'pyproclust.algorithms.hierarchical',
                'pyproclust.algorithms.spectral',
                'pyproclust.algorithms.dbscan',
                'pyproclust.algorithms.kmedoids',
                'pyproclust.algorithms.random',
                'pyproclust.algorithms.gromos',
                'pyproclust.htmlreport',
                'pyproclust.tools',
                'pyproclust.clustering',
                'pyproclust.clustering.analysis',
                'pyproclust.clustering.comparison',
                'pyproclust.clustering.comparison.distrprob',
                'pyproclust.clustering.selection',
                'pyproclust.clustering.metrics',
                'pyproclust.clustering.metrics.cython',
                'pyproclust.clustering.filtering',
                'pyproclust.protocol',
                'pyproclust.protocol.refinement',
                'pyproclust.protocol.exploration',
                'pyproclust.driver',
                'pyproclust.driver.handlers',
                'pyproclust.driver.handlers.matrix',
                'pyproclust.driver.observer',
                'pyproclust.driver.scheduling',
                'pyproclust.driver.compressor'
      ],
      ext_modules=[
                   Extension('pyproclust.clustering.metrics.cython.normNCut',[
                                'pyproclust/clustering/metrics/cython/normNCut.c'                         
                   ]),
                   Extension('pyproclust.clustering.metrics.cython.boundedCohesion', [
                                'pyproclust/clustering/metrics/cython/boundedCohesion.c'                                                 
                   ]),
                   Extension('pyproclust.clustering.metrics.cython.silhouette',[
                                'pyproclust/clustering/metrics/cython/silhouette.c'                                       
                   ]),
                   Extension('pyproclust.clustering.metrics.cython.meanMinimumDistance', [
                                'pyproclust/clustering/metrics/cython/meanMinimumDistance.c'                                                        
                   ])
      ],
      
     )
