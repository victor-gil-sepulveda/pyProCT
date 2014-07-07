'''
Created on 25/02/2013

@author: victor
'''

if __name__ == '__main__': # Comp. with sphynx
    from distutils.core import setup, Extension
    import numpy
    import distutils.sysconfig
    import os

    def read(fname):
        return open(os.path.join(os.path.dirname(__file__), fname)).read()

    setup(
          name='pyProCT',
          version='1.2.0',
          description='pyProCT is an open source cluster analysis software especially adapted for jobs related with structural proteomics',
          author='Victor Alejandro Gil Sepulveda',
          author_email='victor.gil.sepulveda@gmail.com',
          url='https://github.com/victor-gil-sepulveda/pyProCT',
          license = 'LICENSE.txt',
          long_description = read('README.rst'),
          packages=[
                    'pyproct',
                    'pyproct.algorithms',
                    'pyproct.algorithms.hierarchical',
                    'pyproct.algorithms.spectral',
                    'pyproct.algorithms.spectral.cython',
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
                    'pyproct.driver.postprocessing'
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
                       ],extra_compile_args=["-O3","-ffast-math"]),
                       Extension("pyproct.algorithms.spectral.cython.spectralTools", [
                                    'pyproct/algorithms/spectral/cython/spectralTools.c'
                       ],extra_compile_args=["-O3","-ffast-math"])
          ],

          install_requires=[
            "numpy>=1.6.1",
            "PIL>=1.1.6",
            "scipy>=0.9.0",
            "matplotlib>=1.1.1rc",
            "ProDy>=1.4.2",
            "fastcluster>=1.1.6",
            "scikit-learn>=0.12",
            "pyScheduler>=0.1.0",
            "pyRMSD>=4.0.0",
            "mpi4py>=1.3"
          ]
    )
