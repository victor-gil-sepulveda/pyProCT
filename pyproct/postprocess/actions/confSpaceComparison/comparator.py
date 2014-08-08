"""
Created on 12/06/2012

@author: victor
"""
import numpy
from pyproct.postprocess.actions.confSpaceComparison.overlapCalculator import OverlapCalculator
from pyproct.postprocess.actions.confSpaceComparison.tools import mergeSeparatedClusters,\
    getAllElements, calculate_mean_center_differences, calculate_distance_stats

class Separator(object):

    def __init__(self):
        """
        """
        pass

    @classmethod
    def classify(cls, decomposed_clusters):
        """
        Classifies decomposed clusters in 2 categories:
        - pure clusters: Are those that only contain elements from one trajectory.
        - mixed clusters: Are those that contain elements from more than one source trajectory.
        @param decomposed_clusters: A list of decomposed clusters.

        @return: A dictionary with two entries (one per category). Each entry holds a list with the
        decomposed clusters in that category.
        """
        classification = {
                          "pure":{},
                          "mixed":{}
                          }

        for cluster_id in decomposed_clusters:
            decomposed_cluster = decomposed_clusters[cluster_id]

            if len(decomposed_cluster.keys())>1:
                classification["mixed"][cluster_id] = decomposed_cluster
            else:
                classification["pure"][cluster_id] = decomposed_cluster
        return classification

    @classmethod
    def decompose(cls, clusters, traj_ranges):
        """
        Converts the clusters given as input into "decomposed clusters",i.e. a new data structure that divides
        the elements of the cluster according to its source trajectory.
        Example:
        cluster 1 has the following elements:
        [1, 3, 4, 7, 8, 12]
        and traj_ranges = {"traj_A":(0,6),"traj_B":(7,15)
        The resulting "decomposed cluster" will be:
        dc = {"traj_A":[1, 3, 4], "traj_B":[7, 8, 12]}

        @param clusters: A list containing cluster objects.
        @param traj_ranges: A dictionary that contains the starting and ending frame of each trajectory (indexed by
        a trajectory id). The numbering is accumulative, so if we have 2 trajectories of 10 models each, first will
        start in 0 and end in 9 and second will start in 10 and end in 19.

        @return: The list of decomposed clusterings.
        """
        set_ranges = {}
        for traj_id in traj_ranges:
            start,end = traj_ranges[traj_id]
            set_ranges[traj_id] = set(range(start,end+1))

        decomposed_clusters = {}
        for c in clusters:
            elements_set = set(c.all_elements)
            decomposed_cluster = {}
            for traj_id in set_ranges:
                intersection = elements_set.intersection(set_ranges[traj_id])
                if len(intersection) > 0:
                    # Then cluster has elements of this trajectory
                    decomposed_cluster[traj_id] = list(intersection)
            decomposed_clusters[c.id] = decomposed_cluster
        return decomposed_clusters

    @classmethod
    def separate(cls, clusters, traj_ranges):
        """
        Decomposes and separates all the clusters in clustering into pure or mixed clusters depending if their
        elements come from one or more source trajectories.
        @param clusters: A list containing cluster objects.
        @param traj_ranges: A dictionary that contains the starting and ending frame of each trajectory (indexed by
        a trajectory id). The numbering is accumulative, so if we have 2 trajectories of 10 models each, first will
        start in 0 and end in 9 and second will start in 10 and end in 19.

        @return: A dictionary with two entries (one per category). Each entry holds a list with the
        decomposed clusters in that category.
        """
        return cls.classify(cls.decompose(clusters, traj_ranges))

class Analyzer(object):
    def __init__(self):
        """
        """
        pass

    @classmethod
    def run(self, decomposed_clusters, matrix):
        """
        Performs a series of analysis to the whole dataset and each of the decomposed clusters.
        @param decomposed_clusters: A dictionary of decomposed clusters containing 2 keys: "pure" and "mixed",
        each with a list containing decomposed clusters of each kind.
        @param matrix: The work distance matrix.

        @return: The analysis dictionary with all the values.
        """
        analysis = {}
        analysis["total_num_elements"] = 0
        analysis["total_num_clusters"] = 0

        self.analyze_clustering(decomposed_clusters, matrix, analysis)

        self.analyze_clusters(decomposed_clusters, matrix, analysis)

        return analysis

    @classmethod
    def analyze_clustering(cls, separated_decomposed_clusters, distance_matrix, analysis):
        analysis["total_num_clusters"] = 0
        analysis["total_num_elements"] = 0
        analysis["overlap"] = OverlapCalculator.calculate_global_overlap(mergeSeparatedClusters(separated_decomposed_clusters), distance_matrix, 2, 1)
        for cluster_type in separated_decomposed_clusters:
            analysis["num_" + cluster_type] = len(separated_decomposed_clusters[cluster_type])
            analysis["total_num_clusters"] += analysis["num_" + cluster_type]
            analysis["num_" + cluster_type + "_elements"] = numpy.sum([len(getAllElements(separated_decomposed_clusters[cluster_type][dc_id])) for dc_id in separated_decomposed_clusters[cluster_type]])
            analysis["total_num_elements"] += analysis["num_" + cluster_type + "_elements"]
        return cluster_type

    @classmethod
    def analyze_clusters(cls, separated_decomposed_clusters, distance_matrix, analysis):
        for cluster_type in separated_decomposed_clusters:
            for cluster_id in separated_decomposed_clusters[cluster_type]:
                decomposed_cluster = separated_decomposed_clusters[cluster_type][cluster_id]
                analysis[cluster_id] = {"components":decomposed_cluster.keys(),"global":{}}
                analysis[cluster_id]["global"]["mean"], analysis[cluster_id]["global"]["std"], analysis[cluster_id]["global"]["max"] = calculate_distance_stats(getAllElements(decomposed_cluster), distance_matrix)
                analysis[cluster_id]["global"]["num_elements"] = len(getAllElements(decomposed_cluster))

                for traj_id in decomposed_cluster:
                    analysis[cluster_id]["global"][traj_id] = {}
                    analysis[cluster_id]["global"][traj_id]["mean"], analysis[cluster_id]["global"][traj_id]["std"], analysis[cluster_id]["global"][traj_id]["max"] = calculate_distance_stats(decomposed_cluster[traj_id], distance_matrix)
                    analysis[cluster_id]["global"][traj_id]["num_elements"] = len(decomposed_cluster[traj_id])

                if cluster_type == "mixed":
                    analysis[cluster_id]["centers_mean_diff"] = calculate_mean_center_differences(decomposed_cluster, distance_matrix)
                    analysis[cluster_id]["global"]["overlap"] = OverlapCalculator.calculate_cluster_overlap(2, decomposed_cluster, distance_matrix)

