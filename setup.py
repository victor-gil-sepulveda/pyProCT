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
          version='1.7.1',
          description='pyProCT is an open source cluster analysis software especially adapted for jobs related with structural proteomics',
          author='Victor Alejandro Gil Sepulveda',
          author_email='victor.gil.sepulveda@gmail.com',
          url='https://github.com/victor-gil-sepulveda/pyProCT',
          license = 'LICENSE.txt',
          long_description = read('README.rst'),
          packages=[
                    'pyproct',
                    'pyproct.clustering',
                    'pyproct.clustering.algorithms',
                    'pyproct.clustering.algorithms.dbscan',
                    'pyproct.clustering.algorithms.dbscan.cython',
                    'pyproct.clustering.algorithms.gromos',
                    'pyproct.clustering.algorithms.hierarchical',
                    'pyproct.clustering.algorithms.kmedoids',
                    'pyproct.clustering.algorithms.random',
                    'pyproct.clustering.algorithms.spectral',
                    'pyproct.clustering.algorithms.spectral.cython',
                    'pyproct.clustering.evaluation',
                    'pyproct.clustering.evaluation.analysis',
                    'pyproct.clustering.evaluation.metrics',
                    'pyproct.clustering.evaluation.metrics.cython',
                    'pyproct.clustering.evaluation.metrics.cython.graph',
                    'pyproct.clustering.filtering',
                    'pyproct.clustering.protocol',
                    'pyproct.clustering.protocol.exploration',
                    'pyproct.clustering.protocol.refinement',
                    'pyproct.clustering.selection',
                    'pyproct.data',
                    'pyproct.data.handler',
                    'pyproct.data.handler.featurearray',
                    'pyproct.data.handler.protein',
                    'pyproct.data.matrix',
                    'pyproct.data.matrix.featurearray',
                    'pyproct.data.matrix.combination',
                    'pyproct.data.matrix.protein',
                    'pyproct.data.matrix.protein.cases',
                    'pyproct.driver',
                    'pyproct.driver.observer',
                    'pyproct.driver.results',
                    'pyproct.driver.scheduling',
                    'pyproct.driver.time',
                    'pyproct.driver.workspace',
                    'pyproct.postprocess',
                    'pyproct.postprocess.actions',
                    'pyproct.postprocess.actions.confSpaceComparison',
                    'pyproct.tools'
                    
          ],

          include_dirs = [numpy.get_include(),
                          distutils.sysconfig.get_python_inc()],
          ext_modules=[
                       # Graph metrics
                       Extension('pyproct.clustering.evaluation.metrics.cython.graph.nCut',[
                                    'pyproct/clustering/evaluation/metrics/cython/graph/nCut.c'
                       ], extra_compile_args=["-O3","-ffast-math"]),
                       Extension('pyproct.clustering.evaluation.metrics.cython.graph.ratioCut',[
                                    'pyproct/clustering/evaluation/metrics/cython/graph/ratioCut.c'
                       ], extra_compile_args=["-O3","-ffast-math"]),
                       Extension('pyproct.clustering.evaluation.metrics.cython.graph.minMaxCut',[
                                    'pyproct/clustering/evaluation/metrics/cython/graph/minMaxCut.c'
                       ], extra_compile_args=["-O3","-ffast-math"]),
                       Extension('pyproct.clustering.evaluation.metrics.cython.graph.tools',[
                                    'pyproct/clustering/evaluation/metrics/cython/graph/tools.c'
                       ], extra_compile_args=["-O3","-ffast-math"]),
                       
                       # Other metrics
                       Extension('pyproct.clustering.evaluation.metrics.cython.cohesion', [
                                    'pyproct/clustering/evaluation/metrics/cython/cohesion.c'
                       ], extra_compile_args=["-O3","-ffast-math"]),
                       Extension('pyproct.clustering.evaluation.metrics.cython.silhouette',[
                                    'pyproct/clustering/evaluation/metrics/cython/silhouette.c'
                       ], extra_compile_args=["-O3","-ffast-math"]),
                       
                       # Algorithm tools
                       Extension("pyproct.clustering.algorithms.dbscan.cython.cythonDbscanTools", [
                                    'pyproct/clustering/algorithms/dbscan/cython/cythonDbscanTools.c'
                       ],extra_compile_args=["-O3","-ffast-math"]),
                       Extension("pyproct.clustering.algorithms.spectral.cython.spectralTools", [
                                    'pyproct/clustering/algorithms/spectral/cython/spectralTools.c'
                       ],extra_compile_args=["-O3","-ffast-math"])
          ],

          install_requires=[
            "pyRMSD>=4.0.0",
            "pyScheduler>=0.1.0",
            "fastcluster>=1.1.6",
            "ProDy>=1.4.2",
            "numpy>=1.6.1",
            "scipy>=0.9.0",
            "scikit-learn>=0.12",
            "Pillow>=2.6.2",
            "matplotlib>=1.1.1rc",
            "mpi4py>=1.3"
          ]
    )
