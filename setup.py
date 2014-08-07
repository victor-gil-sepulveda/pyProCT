"""
Created on 25/02/2013

@author: victor
"""

if __name__ == '__main__': # Compatibility with sphynx
    from distutils.core import setup, Extension
    import numpy
    import distutils.sysconfig
    import os

    def read(fname):
        return open(os.path.join(os.path.dirname(__file__), fname)).read()

    setup(
          name='pyProCT',
          version='1.4.0',
          description='pyProCT is an open source cluster analysis software especially adapted for jobs related with structural proteomics',
          author='Victor Alejandro Gil Sepulveda',
          author_email='victor.gil.sepulveda@gmail.com',
          url='https://github.com/victor-gil-sepulveda/pyProCT',
          license = 'LICENSE.txt',
          long_description = read('README.rst'),
          packages=[
                    'pyproct',
                    'pyproct.data',
                    'pyproct.data.proteins',
                    'pyproct.data.proteins.matrix',
                    'pyproct.tools',
                    'pyproct.clustering',
                    'pyproct.clustering.algorithms',
                    'pyproct.clustering.algorithms.hierarchical',
                    'pyproct.clustering.algorithms.spectral',
                    'pyproct.clustering.algorithms.spectral.cython',
                    'pyproct.clustering.algorithms.dbscan',
                    'pyproct.clustering.algorithms.dbscan.cython',
                    'pyproct.clustering.algorithms.kmedoids',
                    'pyproct.clustering.algorithms.random',
                    'pyproct.clustering.algorithms.gromos',
                    'pyproct.clustering.evaluation',
                    'pyproct.clustering.evaluation.analysis',
                    'pyproct.clustering.evaluation.metrics',
                    'pyproct.clustering.evaluation.metrics.cython',
                    'pyproct.clustering.selection',
                    'pyproct.clustering.filtering',
                    'pyproct.clustering.protocol',
                    'pyproct.clustering.protocol.refinement',
                    'pyproct.clustering.protocol.exploration',
                    'pyproct.driver',
                    'pyproct.driver.time',
                    'pyproct.driver.observer',
                    'pyproct.driver.scheduling',
                    'pyproct.driver.workspace',
                    'pyproct.driver.results',
                    'pyproct.postprocess',
                    'pyproct.postprocess.actions',
                    'pyproct.postprocess.actions.confSpaceComparison'
                    
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
