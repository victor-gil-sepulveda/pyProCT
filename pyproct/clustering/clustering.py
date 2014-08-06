"""
Created on 20/03/2012

@author: victor
"""
import os
import pickle
import math
from pyproct.tools import scriptTools
from pyproct.clustering.cluster import Cluster

class Clustering(object):
    """
    Class representing the result of a clustering algorithm (a set of clusters).
    """
    def __init__(self, clusters, details="", sort = True):
        """
        Constructor
        """
        self.total_number_of_elements  = 0
        self.clusters = []
        for c in clusters:
            self.add_cluster(c)
        self.sorted = False
        if sort:
            self.sort_clusters_by_size()
        self.details = details

    def add_cluster(self, cluster):
        """
        Adds a cluster to the clusterization.
        """
        self.sorted = False
        # Add cluster_id
        cluster.id = "cluster_"+str(len(self.clusters))
        self.clusters.append(cluster)
        self.total_number_of_elements += cluster.get_size()

    def get_population_percent_of_cluster(self, i):
        """
        Returns the percentage of elements of this cluster over the number of
        all the elements of the dataset used to get the clusters (total_elements).
        """
        if len(self.clusters) > 0 :
            if i < len(self.clusters):
                data_in_this_cluster = self.clusters[i].get_size()
                return data_in_this_cluster * 100. / self.total_number_of_elements
            else:
                print "[WARNING Clustering::get_population_percent_of_cluster] you want the percentage of a non existing cluster."
                return 0
        else:
            print "[WARNING Clustering::get_population_percent_of_cluster] No clusters in this clusterization."
            return 0.

    def get_population_percent_of_n_bigger_clusters(self, n):
        """
        Gets the percentage of population over total population  of the first
        n bigger modules.
        """
        percents = []
        if n > len(self.clusters):
            print "[WARNING get_population_percent_of_n_bigger_clusters] Can't get more percentages than clusters we have. N defaulted to",len(self.clusters)
            n = len(self.clusters)

        # Precondition: clusters are sorted
        self.sort_clusters_by_size()

        for i in range(n):
            percents.append(self.get_population_percent_of_cluster(i))

        return percents

    def sort_clusters_by_size(self, reverse = True):
        """
        Returns a list o sorted clusters by size of the cluster. The biggest cluster goes first
        unless reverse = False.
        """
        if not self.sorted:
            self.clusters.sort(reverse=reverse, key=lambda c: c.get_size())
            self.sorted = True

    def number_of_clusters_to_get_percent(self, percent):
        """
        Returns the number of clusters needed to get the percentage 'percent' over total_elements.
        """
        total_cluster_frame_percent = self.total_number_of_elements* (percent/100.)
        total_frames_visited = 0
        clusters_visited = 0

        # Precondition: clusters are sorted
        self.sort_clusters_by_size()

        for c in self.clusters:
            clusters_visited = clusters_visited + 1
            total_frames_visited = total_frames_visited + c.get_size()
            if total_frames_visited >= total_cluster_frame_percent:
                break

        return clusters_visited

    def cluster_is_inside(self, cluster):
        """
        Returns True if the cluster 'cluster' is currently in the clusters list.
        """
        return self.cluster_index(cluster) != -1

    def cluster_index(self, cluster):
        """
        Returns the index where the cluster 'cluster' is in the list 'self.clusters'
        or -1 if it's not.
        """
        for i in range(len(self.clusters)):
            if cluster == self.clusters[i]:
                return i
        return -1

    def eliminate_cluster(self, cluster):
        """
        Deletes one cluster of the clusters list. The cluster has to be inside the list.
        """
        self.clusters.remove(cluster)
        self.total_number_of_elements -= cluster.get_size()

    def eliminate_noise(self, cluster_size):
        """
        Deletes all the clusters with a number of elements under cluster_size.
        """
        to_eliminate = []
        for c in self.clusters:
            if c.get_size()<cluster_size:
                to_eliminate.append(c)

        for c in to_eliminate:
            self.eliminate_cluster(c)

    def gen_class_list(self, number_of_elements = None, starts_with = 0):
        """
        Generates a class list from a clustering.
        """
        if number_of_elements is not None:
            class_list = [-1]*number_of_elements
        else:
            class_list = [-1]*(max(self.get_all_clustered_elements())+1)
        j = starts_with
        for c in self.clusters:
            for n in c.all_elements:
                class_list[n] = j
            j = j + 1  # one class for each cluster
        return class_list

    def get_all_clustered_elements(self):
        """
        Returns a list with the elements in all its clusters.
        """
        all_clustered_elements = []
        for c in self.clusters:
            all_clustered_elements.extend(c.all_elements)
        return all_clustered_elements

    def save_to_disk(self, path_with_file_name):
        """
        Saves itself as a binary file.

        @param filename: complete path with name of the file
        """
        not_repeated_path_with_file_name = scriptTools.get_not_repeated_file_name(path_with_file_name)
        file_handler = open(not_repeated_path_with_file_name,'w')
        pickle.dump(self,file_handler)
        file_handler.close()

    @classmethod
    def load_from_disk(cls, path_with_file_name):
        """
        Loads a clustering which was stored into disk.

        @param path_with_file_name: complete path with name of the file (the same used for saving).

        @return: The clustering stored in this file.
        """
        file_handler = open(path_with_file_name,'r')
        clustering = pickle.load(file_handler)
        file_handler.close()
        return clustering

    @classmethod
    def load_all_from_directory(cls, directory):
        """
        Loads all clusterings residing in a directory.

        @param directory: The directory path.

        @return: The list of loaded clusterings (with their filenames) inside this directory.
        """
        clusterings = []
        clustering_files = []

        files = os.listdir(directory)
        for filename in files:
            if ".bin" in filename:
                clustering_files.append(directory+"/"+filename)
        clustering_files.sort()

        for a_file in clustering_files:
            clusterings.append((cls.load_from_disk(a_file), a_file))

        return clusterings

    @classmethod
    def classify(cls,tags,clusterings):
        """
        Classifies a collection of clusterings using the classes in 'tags', counting the occurrences for each class.
        A clustering belongs to a class if in its 'details' string the class-tag appears.

        @param tags: An array with the class-tags used to classify.
        @param clusterings: The list of clusterings to classify.

        @return: A counter of the class occurrences.
        """
        counter = {}
        for t in tags:
            counter[t] = 0
        for clustering in clusterings:
            for t in tags:
                if t in clustering.details:
                    counter[t] += 1
        return counter

    def get_medoids(self, distance_matrix):
        """
        Returns a list containing the medoids of all clusters.

        @param distance_matrix: Is the distance matrix used to calculate the clustering.

        @return: The list of medoids.
        """
        medoids= []
        for c in self.clusters:
            medoids.append(c.calculate_medoid(distance_matrix))
        return medoids

    def get_proportional_size_representatives(self, number_of_structures, distance_matrix):
        """
        Returns a list with the medoids and a random sample of cluster elements so that the number of elements
        of each cluster in the list is proportional to the size of each cluster.

        @param number_of_structures: Is the final number of elements we want in our list.

        @param distance_matrix: Is the distance matrix used to calculate the clustering.

        @return: The list of representative elements.
        """
        representatives = []
        for c in self.clusters:
            number_of_elements = max(0,int((math.ceil(c.get_size() / float(self.total_number_of_elements)*number_of_structures)))-1) # minus the medoid
            representatives.extend( c.get_random_sample(number_of_elements))
        representatives.extend(self.get_medoids(distance_matrix))
        return representatives

    def to_dic(self):
        """
        Converts this clustering into a dictionary (to be used with json serializers).
        """
        clustering_dic = {}
        clustering_dic["total_number_of_elements"] = self.total_number_of_elements
        clustering_dic["number_of_clusters"] = len(self.clusters)
        clustering_dic["clusters"] = []
        for cluster in self.clusters:
            clustering_dic["clusters"].append(cluster.to_dic())
        return clustering_dic

    @classmethod
    def from_dic(cls, clustering_dic):
        """
        Creates a clustering from a clustering dictionary describing it (as reverse operation of
        'to_dic').
        """
        clusters_dic = clustering_dic["clusters"];
        clusters = []
        for cluster_dic in clusters_dic:
            clusters.append(Cluster.from_dic(cluster_dic))

        return Clustering(clusters)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            my_dic_clusters_elems = [x["elements"] for x in self.to_dic()["clusters"]]
            other_dic_clusters_elems = [x["elements"] for x in other.to_dic()["clusters"]]
            if len(my_dic_clusters_elems) != len(other_dic_clusters_elems):
                return False
            else:
                for cluster_elems in my_dic_clusters_elems:
                    if not cluster_elems in other_dic_clusters_elems:
                        return False
                return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
