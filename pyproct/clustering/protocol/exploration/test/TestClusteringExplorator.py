"""
Created on 07/02/2013

@author: victor
"""
import unittest
from pyproct.clustering.protocol.exploration.clusteringExplorer import ClusteringExplorer


class TestClusteringExplorator(unittest.TestCase):

    def test_generate_clustering_info(self):
        explorator = ClusteringExplorer({"clustering":{"evaluation":{}}},None, None, None, None)
        run_parameters = []
        for i in range(10):
            run_parameters.append({"lol":i})

        expected = {
                    'clustering_0000': {'clustering': None, 'type': 'my_algorithm', 'parameters': {"lol":0}},
                    'clustering_0001': {'clustering': None, 'type': 'my_algorithm', 'parameters': {"lol":1}},
                    'clustering_0002': {'clustering': None, 'type': 'my_algorithm', 'parameters': {"lol":2}},
                    'clustering_0003': {'clustering': None, 'type': 'my_algorithm', 'parameters': {"lol":3}},
                    'clustering_0004': {'clustering': None, 'type': 'my_algorithm', 'parameters': {"lol":4}},
                    'clustering_0005': {'clustering': None, 'type': 'my_algorithm', 'parameters': {"lol":5}},
                    'clustering_0006': {'clustering': None, 'type': 'my_algorithm', 'parameters': {"lol":6}},
                    'clustering_0007': {'clustering': None, 'type': 'my_algorithm', 'parameters': {"lol":7}},
                    'clustering_0008': {'clustering': None, 'type': 'my_algorithm', 'parameters': {"lol":8}},
                    'clustering_0009': {'clustering': None, 'type': 'my_algorithm', 'parameters': {"lol":9}},
                    }
        clustering_info = explorator.generate_clustering_info("my_algorithm", run_parameters)
        self.assertEqual(clustering_info, expected)

        expected = {
                    'clustering_0010': {'clustering': 0, 'type': 'my_algorithm', 'parameters': {'lol': 0}},
                    'clustering_0011': {'clustering': 1, 'type': 'my_algorithm', 'parameters': {'lol': 1}},
                    'clustering_0012': {'clustering': 2, 'type': 'my_algorithm', 'parameters': {'lol': 2}},
                    'clustering_0013': {'clustering': 3, 'type': 'my_algorithm', 'parameters': {'lol': 3}},
                    'clustering_0014': {'clustering': 4, 'type': 'my_algorithm', 'parameters': {'lol': 4}},
                    'clustering_0015': {'clustering': 5, 'type': 'my_algorithm', 'parameters': {'lol': 5}},
                    'clustering_0016': {'clustering': 6, 'type': 'my_algorithm', 'parameters': {'lol': 6}},
                    'clustering_0017': {'clustering': 7, 'type': 'my_algorithm', 'parameters': {'lol': 7}},
                    'clustering_0018': {'clustering': 8, 'type': 'my_algorithm', 'parameters': {'lol': 8}},
                    'clustering_0019': {'clustering': 9, 'type': 'my_algorithm', 'parameters': {'lol': 9}}
                    }

        clustering_info = explorator.generate_clustering_info("my_algorithm", run_parameters, range(10))
        self.assertEqual(clustering_info, expected)

    def test_build_algorithm(self):
        self.fail("TODO")

    def test_run(self):
        self.fail("TODO")

    def test_schedule_algorithm(self):
        self.fail("TODO")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()