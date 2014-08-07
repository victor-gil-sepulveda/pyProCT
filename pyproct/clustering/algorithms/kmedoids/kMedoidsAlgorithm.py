"""
Created on 29/05/2012

@author: victor
"""

import random
from pyproct.clustering.algorithms.gromos.gromosAlgorithm import GromosAlgorithm
from pyproct.clustering.cluster import gen_clusters_from_class_list
from pyproct.clustering.clustering import Clustering



class KMedoidsAlgorithm(object):
    """
    K-Means like algorithm with medoids (so that the cluster prototype is always a member of the
    dataset). It has 3 different initial seeding implementations.
    """

    MAX_ITERATIONS = 500

    @staticmethod
    def seeding_types():
        """
        Returns the possible initial seeding types.
            - GROMOS: Uses GROMOS to get the initial medoids. Useful if clusters are well separated.
            - EQUIDISTANT: Divides the dataset in k consecutive parts and uses their central element as
                medoid. Useful if we suspect that sequence order and geometrical likeness are correlated
                (for instance with MD sequences).
            - RANDOM: Uses a random choice of elements from the dataset as initial medoids.
        """
        return  ["GROMOS","RANDOM", "EQUIDISTANT"]

    def __init__(self, condensed_matrix, **kwargs):
        """
        Constructor. Sets up some class variables.

        @param condensed_matrix: The dataset's distance matrix.
        @param rand_seed: A number to set the random seed or None if we don't want to set it. Useful
        to be able to reproduce results.
        """
        try:
            rand_seed = kwargs["rand_seed"]
            random.seed(rand_seed)
        except KeyError:
            random.seed()

        self.condensed_matrix = condensed_matrix
        self.total_elements = condensed_matrix.row_length
        self.class_list = [0]*condensed_matrix.row_length
        self.gromos_clusters_bookkeeping = {}

    def perform_clustering(self, kwargs):
        """
        Does the actual clustering.
        @param kwargs: Dictionary with this mandatory keys:
            - 'k': Number of clusters to generate.
            - 'seeding_type': One of the initial medoid selectors available (@see seeding_types() ).
                If seeding type is 'GROMOS', 'seeding_max_cutoff' must be also defined, containing the
                cutoff that the GROMOS Algorithm will use. Default is EQUIDISTANT
        """
        self.k = kwargs["k"]

        self.seeding_type =  kwargs["seeding_type"] if "seeding_type" in kwargs else  "EQUIDISTANT"

        if self.seeding_type == 'GROMOS':
            self.seeding_max_cutoff = kwargs["seeding_max_cutoff"]
        else:
            self.seeding_max_cutoff = -1.0

        # Getting first medoids
        current_medoids = self.seeding(self.k, self.seeding_max_cutoff, self.seeding_type)
#         print "Initial medoids:", current_medoids

        current_iteration = 0
        last_medoids = []
        while not self.convergence(current_medoids,last_medoids) and current_iteration < KMedoidsAlgorithm.MAX_ITERATIONS:
#             print "Iteration"
            # Update clusters
            self.cluster_update(current_medoids, self.condensed_matrix)
            # Copy current to last (by reference)
            last_medoids = current_medoids
            # Update medoids
            current_medoids = self.update_medoids()
#             print "Current medoids:", current_medoids
            current_iteration = current_iteration + 1

        algorithm_details = "K-Medoids algorithm with k ="+str(int(self.k))+" and %s initial seeding"%self.seeding_to_str()
        clusters = gen_clusters_from_class_list(self.class_list)
        return Clustering(clusters,details = algorithm_details)

    def seeding_to_str(self):
        """
        Convenience function to build the details string.

        @return: a string representation of the initial seeding type.
        """
        seeding_string = self.seeding_type
        if self.seeding_type == 'GROMOS':
            seeding_string = "(GROMOS,"+str(self.seeding_max_cutoff)+")"
        return seeding_string


    def convergence(self, current, last):
        """
        Checks if the medoid list has changed. If not, then the algorithm has converged.

        @param current: The current medoid list in this step.
        @param last: The medoid list of the previous step.

        @return: True if the lists are equal, False otherwise.
        """
#        print "Checking convergence"
        if last == []:
            return False

        current.sort()
        last.sort()
        for i in range(len(current)):
            if current[i] != last[i]:
                return False
        return True

    def get_closer_medoid(self, element, medoids, condensed_matrix):
        """
        Says which medoid is closer to an element.

        @param element: The element of the dataset we want to test.
        @param condensed_matrix: The dataset's distance matrix.

        @return: The closer medoid to this element.
        """
        distances = []
        for j in range(len(medoids)):
            distances.append((condensed_matrix[element,medoids[j]],medoids[j]))
        (winner_dist, winner_medoid) = min(distances)
        return winner_medoid

    def gen_medoid_to_cluster_id_map(self, medoids):
        """
        It generates a map defining the relationship between cluster classes and medoids.
        For instance if k = 3 and our medoid list is [10, 23, 50], class 0 corresponds to medoid 10, so
        cluster_id_map[10] = class 0

        @param medoids: The medoid list.

        @return: A medoid to class dictionary.
        """
        cluster_id_map = {}
        for i in range(len(medoids)):
            cluster_id_map[medoids[i]] = i
        return cluster_id_map

    def cluster_update(self, medoids, condensed_matrix):
        """
        Assigns each element to a cluster (updating its class), depending on which medoid is closer.

        @param medoids: The medoid list.
        @param condensed_matrix: The dataset's distance matrix.

        """
        cluster_id_map = self.gen_medoid_to_cluster_id_map(medoids)
#         print "cluster medoid -> class map:", cluster_id_map
        for i in range(self.total_elements):
            self.class_list[i] = cluster_id_map[self.get_closer_medoid(i,medoids,condensed_matrix)]

    def update_medoids(self):
        """
        Regenerates the medoids list once the new clusters have been generated.

        @return: A new medoid list.
        """
        clusters = gen_clusters_from_class_list(self.class_list)
        medoids = []
        for c in clusters:
            medoids.append(c.calculate_medoid(self.condensed_matrix))
        return medoids


    def seeding(self, k, seeding_max_cutoff, seeding_type):
        """
        Does the initial seeding using the stated seeding process.

        @param k: The number of medoids we need.
        @param seeding_max_cutoff: The initial seeding cutoff to use with the GROMOS algorithm.
        @param seeding_type: The seeding process type (@see seeding_types() )

        @return: The initial medoids.
        """
        if not seeding_type in self.seeding_types():
            print "[ERROR::SpectralClusteringAlgorithm] Seeding type " ,seeding_type, "is not a correct type. Use one of these instead: ", self.seeding_types()
            exit()

        if seeding_type == "RANDOM":
            return self.random_seeding(k)

        elif seeding_type == "GROMOS":
            return self.gromos_seeding(k, seeding_max_cutoff)

        elif seeding_type == "EQUIDISTANT":
            return self.equidistant_seeding(k, self.condensed_matrix.row_length)

    def random_seeding(self,k):
        """
        Returns k random medoids.

        @param k: The number of medoids we need.

        @return: The medoid list.
        """
        random_medoids = random.sample(range(self.condensed_matrix.row_length),k)
        return random_medoids

    def equidistant_seeding(self, k, number_of_elements):
        """
        @see seeding_types

        @param k: The number of medoids we need.
        @param number_of_elements: Number of elements of the dataset.

        @return: The medoid list.
        """
        step = number_of_elements / k
        medoids = []
        for i in range(k):
            medoids.append((step/2) + (i*step))

        return medoids

    def gromos_seeding(self, k, seeding_max_cutoff):
        """
        The first k medoids are selected running the gromos algorithm with an initial cutoff defined in 'seeding_max_cutoff'
        until it gets at least k clusters. If it can't, it will use the random seeding strategy.

        @param k: The number of medoids we need(and thus the number of clusters we want GROMOS to produce)
        @param seeding_max_cutoff: The initial seeding cutoff to use with the gromos algorithm.

        @return: The medoid list.
        """
#         print "Seeding process"
        current_cutoff = seeding_max_cutoff
        clusters = []
        while len(clusters) < k and current_cutoff > 0:
            """
            EL BOOKEEPING SOLO TIENE SENTIDO EN MULTITHREADING; PERO NO
            CON PROCESOS SEPARADOS :(
            if(current_cutoff in self.gromos_clusters_bookkeeping):
                if
            else:
            """
            print "Trying gromos with cutoff =", current_cutoff,"for seeding"
            gromos_algorithm = GromosAlgorithm(self.condensed_matrix)
            clusters = gromos_algorithm.perform_clustering({"cutoff":current_cutoff, "max_clusters":k}).clusters
#             self.gromos_clusters_bookkeeping[current_cutoff] = len(clusters)
            current_cutoff -= 0.1

        # If it was impossible, do a random seeding
        if(current_cutoff<=0.):
            print "Returning a random sampling."
            return self.random_seeding(k)
        else:
            print "Returning the medoids."
            medoids = []
            for c in clusters[0:k]:
                medoids.append(c.prototype)
            return medoids

