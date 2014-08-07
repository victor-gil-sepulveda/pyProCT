"""
Created on 7/8/2014

@author: victor
"""

if __name__ == '__main__': # Compatibility with sphynx
    from setuptools import setup

    setup(
          name='plugin_test_2',
          packages=[
                    'pack2'
          ],
          entry_points = """
                [pyproct.test.plugin]
                testplugin1 = pack2.packmodule:handle_matrix_test_class"""
    )
