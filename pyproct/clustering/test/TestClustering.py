"""
Created on 20/03/2012

@author: victor
"""
import unittest
from pyproct.clustering.cluster import Cluster
from pyproct.clustering.clustering import Clustering
import pyproct.tools.test.data as test_data
import cStringIO
import os
import numpy

class ClusterMock():
    def __init__(self, elements):
        self.elements = elements

    def calculate_medoid(self, param):
        return self.elements[0]

    def get_random_sample(self, n):
        return self.elements[0:n]

    def get_size(self):
        return len(self.elements)

class TestClustering(unittest.TestCase):

    def test_get_percent_population_of_cluster(self):
        clusters =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization = Clustering(clusters)

        total = 0
        for i in range(4):
            total = total + clusterization.get_population_percent_of_cluster(i)
        self.assertAlmostEqual(total, 100., 2)

    def test_get_percent_of_n_clusters(self):
        clusters =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization = Clustering(clusters)

        percents = clusterization.get_population_percent_of_n_bigger_clusters(3)
        expected_percents = [41.1764705882,29.4117647059,23.5294117647]
        for i in range(3):
            self.assertAlmostEqual(percents[i], expected_percents[i], 1)

    def test_number_of_clusters_needed_to_get_this_percent_of_elems(self):
        clusters =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )

        clusterization = Clustering(clusters)

        self.assertEqual(clusterization.number_of_clusters_to_get_percent( 71),3)
        self.assertEqual(clusterization.number_of_clusters_to_get_percent( 70),2)
        self.assertEqual(clusterization.number_of_clusters_to_get_percent( 40),1)
        self.assertEqual(clusterization.number_of_clusters_to_get_percent( 42),2)
        self.assertEqual(clusterization.number_of_clusters_to_get_percent( 100),4)

    def test_creation(self):
        # The inner list is a copy but shares clusters
        clusters =(
                      Cluster(16,[16]),
                      Cluster(4,[4,5,6,7,8]),
                      Cluster(0,[0,1,2,3]),
                      Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization = Clustering(clusters)
        clusters[1].prototype = -20
        self.assertEqual(clusters[1].prototype, clusterization.clusters[1].prototype)

    def test_cluster_is_inside(self):
        clusters =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        not_in_cluster= Cluster(17,[17,16])
        in_cluster = Cluster(0,[0,1,2,3])
        clusterization = Clustering(clusters)
        self.assertEqual(clusterization.cluster_index(not_in_cluster),-1)
        self.assertEqual(clusterization.cluster_index(in_cluster),2)
        self.assertEqual(clusterization.cluster_is_inside(not_in_cluster),False)
        self.assertEqual(clusterization.cluster_is_inside(in_cluster),True)

    def test_remove_cluster(self):
        clusters =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization = Clustering(clusters)
        c = Cluster(0,[0,1,2,3])
        clusterization.eliminate_cluster(c)
        self.assertEqual(len(clusterization.clusters), 3)

    def test_remove_noise(self):
        clusters =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization = Clustering(clusters)
        clusterization.eliminate_noise(5)
        self.assertEqual(len(clusterization.clusters), 2)

    def test_gen_class_list(self):
        clusters =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clustering = Clustering(clusters)
        class_list = clustering.gen_class_list()
        expected_class_list = [2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 3]
        self.assertItemsEqual(class_list, expected_class_list)

        clusters =(
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization = Clustering(clusters)
        class_list = clusterization.gen_class_list()
        expected_class_list = [1, 1, 1, 1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0]
        self.assertItemsEqual(class_list, expected_class_list)

    def test_get_all_clustered_elements(self):
        clusters =(
                  Cluster(16,[16]),
                  Cluster(4,[4,5,6,7,8]),
                  Cluster(0,[0,1,2,3]),
                  Cluster(9,[9,10,11,12,13,14,15])
                  )
        clusterization = Clustering(clusters)
        self.assertItemsEqual(sorted( clusterization.get_all_clustered_elements()), range(17))

    def test_load_and_save_to_disk(self):
        clusters =(Cluster(16,[16]),
                   Cluster(4,[4,5,6,7,8]),
                   Cluster(0,[0,1,2,3]),
                   Cluster(9,[9,10,11,12,13,14,15]))

        clustering = Clustering(clusters)
        before_saving_elements = clustering.get_all_clustered_elements()
        clustering.save_to_disk("data/saved_clustering_for_test")
        loaded_clustering = Clustering.load_from_disk("data/saved_clustering_for_test")
        after_saving_elements = loaded_clustering.get_all_clustered_elements()
        self.assertItemsEqual(before_saving_elements, after_saving_elements)
        os.system("rm data/saved_clustering_for_test")

    def test_batch_load(self):
        clusters =((Cluster(16,[16]), "data/training_clustering_1.bin"),
                   (Cluster(4,[4,5,6,7,8]), "data/training_clustering_2.bin"),
                   (Cluster(0,[0,1,2,3]), "data/training_clustering_3.bin"),
                   (Cluster(9,[9,10,11,12,13,14,15]), "data/training_clustering_4.bin"))

        # Creates 4 clusterings of 1 cluster
        filenames = []
        for cluster, filename in clusters:
            Clustering([cluster]).save_to_disk(filename)
            filenames.append(filename)

        # Then loads them and extracts its elements
        elements = []
        for filename in filenames:
            elements.extend(Clustering.load_from_disk(filename).get_all_clustered_elements())

        elements_batch = []
        clusterings_batch = Clustering.load_all_from_directory("data/")
        for clustering, filename in clusterings_batch:
            elements_batch.extend(clustering.get_all_clustered_elements())

        # And cleans the house
        for filename in filenames:
            os.system("rm "+filename)

        numpy.testing.assert_equal(sorted(elements), range(17))
        numpy.testing.assert_equal(sorted(elements_batch), range(17))

    def test_classify(self):
        tags = ["A","B","C"]
        clusterings = [Clustering([], "this is of type A"),Clustering([], "this is of type B"),Clustering([], "this is of type C"),
                       Clustering([], "this is of type B"),Clustering([], "this is of type S"),Clustering([], "this is of type A"),
                       Clustering([], "this is of type A"),Clustering([], "this is of type C"),Clustering([], "this is of type D")]
        counter =  Clustering.classify(tags, clusterings)
        self.assertEqual(counter['A'], 3)
        self.assertEqual(counter['B'], 2)
        self.assertEqual(counter['C'], 2)

    def test_get_medoids(self):
        clusters = [ClusterMock(range(0,10)),ClusterMock(range(10,50)),ClusterMock(range(50,80)),ClusterMock(range(80,200))]
        clustering = Clustering(clusters)
        self.assertItemsEqual(clustering.get_medoids("distance_matrix"),[0, 10, 50, 80])

    def test_get_proportional_size_representatives(self):
        clusters = [ClusterMock(range(0,10)),ClusterMock(range(10,50)),ClusterMock(range(50,80)),ClusterMock(range(80,200))]
        clustering = Clustering(clusters)
        rep =  clustering.get_proportional_size_representatives(30, "distance_matrix" )
        self.assertItemsEqual(rep, [0, 0, 10, 10, 11, 12, 13, 14, 50, 50, 51, 52, 53, 80, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96])

    def test_to_dic(self):
        clustering =Clustering([Cluster(16,[16]),
                                 Cluster(9,[9,10,11,12,13,14,15]),
                                 Cluster(0,[0,1,2,3]),
                                 Cluster(4,[4,5,6,7,8])])
        self.assertDictEqual(clustering.to_dic(),
                             {'clusters': [{'prototype': 9, 'elements': '9:15', 'id': 'cluster_1'},
                                           {'prototype': 4, 'elements': '4:8', 'id': 'cluster_3'},
                                           {'prototype': 0, 'elements': '0:3', 'id': 'cluster_2'},
                                           {'prototype': 16, 'elements': '16', 'id': 'cluster_0'}],
                              'total_number_of_elements': 17,
                              'number_of_clusters': 4})

    def test_equality(self):
        clusteringA =Clustering([Cluster(16,[16]),
                                 Cluster(4,[4,5,6,7,8]),
                                 Cluster(0,[0,1,2,3]),
                                 Cluster(9,[9,10,11,12,13,14,15])])

        clusteringB =Clustering([Cluster(16,[16]),
                                 Cluster(9,[9,10,11,12,13,14,15]),
                                 Cluster(0,[0,1,2,3]),
                                 Cluster(4,[4,5,6,7,8])])

        clusteringC =Clustering([Cluster(13,[13]),
                                 Cluster(9,[9,10,11,12,16,14,15]),
                                 Cluster(0,[0,1]),
                                 Cluster(4,[2,3,4,5,6,7,8])])

        self.assertEqual(clusteringA, clusteringB)
        self.assertNotEqual(clusteringA, clusteringC)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()