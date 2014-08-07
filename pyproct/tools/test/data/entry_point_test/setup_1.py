"""
Created on 7/8/2014

@author: victor
"""

if __name__ == '__main__': # Compatibility with sphynx
    from setuptools import setup

    setup(
          name='plugin_test',
          packages=[
                    'pack1'
          ],
          entry_points = """
                [pyproct.test.plugin]
                testplugin1 = pack1.packmodule:handle_matrix_test_class
                testplugin2 = pack1.packmodule:handle_analysis_test_class"""
    )
