'''
Created on 12/06/2012

@author: victor
'''
import numpy
from pyproct.clustering.cluster import Cluster
from pyproct.clustering.metrics.common import get_distances_of_elements_to

def calculate_mean_center_differences(decomposed_cluster, matrix):
    """
    Given a mixed decomposed cluster, it calculates the mean of all center differences (giving a qualitative
    view of how separated the inner subclusters are).
    @param decomposed_cluster: A MIXED decomposed cluster.
    @param matrix: The condensed distance matrix used.

    @return: The mean of center distances.
    """
    centers = []
    for traj_id in decomposed_cluster:
        cluster = Cluster(None, decomposed_cluster[traj_id])
        centers.append(cluster.calculate_medoid(matrix))

    center_distances = []
    for i in centers[:-1]:
        for j in centers[1:]:
            center_distances.append(matrix[i,j])
    return numpy.mean(center_distances)

def calculate_distance_stats(elements, matrix):
    """
    Calculates the mean, dispersion and radius of all the distances to the central element of a set of
    elements.
    @param elements: The elements we are working with.
    @param matrix: The used condensed matrix.

    @return: Mean, std deviation and radius of all the elements with respect to their central element.
    """
    cluster = Cluster(None, elements)
    medoid = cluster.calculate_medoid(matrix)

    # We also get a 0 distance from the medoid vs itself (it is contained in 'elements')
    distances = get_distances_of_elements_to(medoid, elements, matrix)

    return numpy.mean(distances), numpy.std(distances), numpy.max(distances)

def getTotalNumberOfElements(list_of_cluster_lists):
    """
    """
    total = 0
    for cluster_list in list_of_cluster_lists:
        for c in cluster_list:
            total += c.get_size()
    return total

class Separator(object):

    def __init__(self):
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
                          "pure":[],
                          "mixed":[]
                          }

        for decomposed_cluster in decomposed_clusters:
            if len(decomposed_clusters.keys())>1:
                classification["mixed"].append(decomposed_cluster)
            else:
                classification["pure"].append(decomposed_cluster)
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
            set_ranges[traj_id].append(set(range(start,end)))

        decomposed_clusters = []
        for c in clusters:
            elements_set = set(c.all_elements)
            decomposed_cluster = {}
            for traj_id in set_ranges:
                intersection = elements_set.intersection(set_ranges[traj_id])
                if len(intersection) > 0:
                    # Then cluster has elements of this trajectory
                    decomposed_cluster[traj_id] = intersection
            decomposed_clusters.append(decomposed_cluster)
        return decomposed_clusters

    @classmethod
    def separate(cls, clusters, traj_ranges):
        """
        Decomposes and separates all the clusters in clustering into pure or mixed clusters depending if their
        elements come from one or more source trajectories.
        """
        return cls.classify(cls.decompose(clusters, traj_ranges))

class ClusterStatisticalData(object):
    def __init__(self, cluster, cluster_type, condensed_distance_matrix,A_elements = [], B_elements = []):
        self.cluster_type = cluster_type
        self.number_of_elements = cluster.get_size()
        self.center_difference = 0
        self.mean_dist_to_center, self.dispersion = calculate_cluster_mean_distance_and_dispersion(cluster,cluster.all_elements,condensed_distance_matrix)
        self.max_distance = calculate_max_distance(cluster,condensed_distance_matrix)

        # If mixed ...
        self.A_part_cluster_data = None
        self.B_part_cluster_data = None
        if cluster_type == 'M':
            self.A_part_cluster = Cluster(None,A_elements)
            self.B_part_cluster = Cluster(None,B_elements)
            self.A_part_cluster_data = ClusterStatisticalData(self.A_part_cluster,'P',condensed_distance_matrix)
            self.B_part_cluster_data = ClusterStatisticalData(self.B_part_cluster,'P',condensed_distance_matrix)
            self.center_difference = calculate_centers_difference(cluster_type, cluster, condensed_distance_matrix,A_elements,B_elements)

    def __repr__(self):
        return str(self)

    def __str__(self):
        res = ""
        if self.cluster_type == 'A' or self.cluster_type == 'B' or self.cluster_type == 'P':
            res = "[Cluster Stats. Data -\n\tType: 'Pure'\n\tSize: %d\n\tMean Dist. to centre: %.3f\n\tDispersion: %.3f]\n"%(self.number_of_elements,\
                                                                                                                     self.mean_dist_to_center,\
                                                                                                                     self.dispersion)
        else:
            res = "[Cluster Stats. Data -\n\ttype: 'Mixed'\n\tSize: %d\n\tMean Dist. to centre: %.3f\n\tDispersion: %.3f\n\tDist. between centers: %.3f\n\n\tA elements inside mixed:"%(self.number_of_elements,\
                                                                                                                     self.mean_dist_to_center,\
                                                                                                                     self.dispersion,self.center_difference)\
            +str(self.A_part_cluster_data )+"\n\tB elements inside mixed:"+str(self.B_part_cluster_data )+"]\n"
        return res

class ClusteringStatisticalAnalyzer(object):
    def __init__(self, clustering, traj_A_pdb, traj_B_pdb, condensed_matrix, trajectory_comparison = True):
        self.traj_A_pdb = traj_A_pdb
        self.traj_B_pdb = traj_B_pdb
        self.matrix = condensed_matrix
        self.clustering = clustering
        self.total_number_of_elements = clustering.total_number_of_elements
        self.trajectory_comparison = trajectory_comparison

        if not trajectory_comparison:
            self.pure_A = clustering.clusters
            self.pure_B = []
            self.mixed_clusters_with_elements = []
        else:
            separator = Separator()
            self.pure_A , self.pure_B, self.mixed_clusters_with_elements = separator.separate(clustering, traj_A_pdb, traj_B_pdb)
            self.traj_A_number_of_models = separator.traj_1_numbr_of_models
            self.traj_B_number_of_models = separator.traj_2_numbr_of_models

        self.trajectory_comparison = trajectory_comparison

        self.cluster_statistical_data = []
        self.statistics_dic = {}

        self.mixed_clusters_without_elements = []
        for m in self.mixed_clusters_with_elements:
            self.mixed_clusters_without_elements.append(m[0])

    def per_cluster_analytics(self):
        for cluster in self.pure_A:
            self.cluster_statistical_data.append(ClusterStatisticalData(cluster,'A',self.matrix))

        for cluster in self.pure_B:
            self.cluster_statistical_data.append(ClusterStatisticalData(cluster,'B',self.matrix))

        for cluster, a_elems, b_elems in self.mixed_clusters_with_elements:
            self.cluster_statistical_data.append(ClusterStatisticalData(cluster,'M',self.matrix,a_elems,b_elems))

    def per_clustering_analytics(self):
        if self.trajectory_comparison:
            self.statistics_dic["number_elements_pure_A"] = getTotalNumberOfElements([self.pure_A])
            self.statistics_dic["elems_percent_pure_A"]= self.statistics_dic["number_elements_pure_A"]*100./self.total_number_of_elements
            self.statistics_dic["elems_self_percent_pure_A"]= self.statistics_dic["number_elements_pure_A"]*100./self.traj_A_number_of_models
            self.statistics_dic["number_elements_pure_B"] = getTotalNumberOfElements([self.pure_B])
            self.statistics_dic["elems_percent_pure_B"]= self.statistics_dic["number_elements_pure_B"]*100./self.total_number_of_elements
            self.statistics_dic["elems_self_percent_pure_B"]= self.statistics_dic["number_elements_pure_B"]*100./self.traj_B_number_of_models
            self.statistics_dic["number_elements_mixed"] = getTotalNumberOfElements([self.mixed_clusters_without_elements])
            self.statistics_dic["elems_percent_mixed"]= self.statistics_dic["number_elements_mixed"]*100./self.total_number_of_elements
