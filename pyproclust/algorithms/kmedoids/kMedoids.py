'''
Created on 29/05/2012

@author: victor
'''

import random
from pyproclust.algorithms.gromos.gromosAlgorithm import GromosAlgorithm
from pyproclust.clustering.cluster import gen_clusters_from_class_list
from pyproclust.clustering.clusterization import Clustering



class KMedoids(object):
    """
    K-Means like algorithm with medoids. Seeding is done with GROMOS.
    """
    def __init__(self,condensed_matrix):
        self.condensed_matrix = condensed_matrix
        self.total_elements = condensed_matrix.row_length
        self.class_list = [0]*condensed_matrix.row_length
        self.gromos_clusters_bookkeeping = {}
        
    def perform_clustering(self, kwargs):    
        self.k = kwargs["k"]
        self.seeding_max_cutoff = kwargs["seeding_max_cutoff"]
        
        # Getting first medoids
        current_medoids = self.seeding(self.k,self.seeding_max_cutoff)
        
        last_medoids = []
        while not self.convergence(current_medoids,last_medoids):
            # Update clusters
            self.cluster_update(current_medoids,self.condensed_matrix)
            # Copy current to last
            last_medoids = current_medoids
            # Update medoids
            current_medoids = self.update_medoids()
    
        algorithm_details = "K-Medoids algorithm with k ="+str(int(self.k))+" and seeding cutoff = "+str(self.seeding_max_cutoff)
        clusters = gen_clusters_from_class_list(self.class_list)
        return Clustering(clusters,details = algorithm_details)
    
    def convergence(self,current, last):
        """
        Checks if the medoid lists have changed. If not, then the algorithm has converged. 
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
            
    
    def get_closer_medoid(self,element,medoids,condensed_matrix):
        """
        Says wich medoid is closer to an element.
        """
        distances = []
        for j in range(len(medoids)):
            distances.append((condensed_matrix[element,medoids[j]],medoids[j]))
        (winner_dist, winner_medoid) = min(distances)
        del winner_dist
        return winner_medoid
    
    def gen_medoid_to_cluster_id_map(self,medoids):
        cluster_id_map = {}
        for i in range(len(medoids)):
            cluster_id_map[medoids[i]] = i
        return cluster_id_map 
    
    def cluster_update(self, medoids, condensed_matrix):
        """
        Assigns each element to a cluster, depending on which medoid is closer.
        """
        cluster_id_map = self.gen_medoid_to_cluster_id_map(medoids)
        for i in range(self.total_elements):
            self.class_list[i] = cluster_id_map[self.get_closer_medoid(i,medoids,condensed_matrix)]
    
    def update_medoids(self):
        """
        Regenerates the medoids list once the new clusters have been done.
        """
        clusters = gen_clusters_from_class_list(self.class_list)
        medoids = []
        for c in clusters:
            medoids.append(c.calculate_medoid(self.condensed_matrix))
        return medoids
    
    def random_seeding(self,k):
        """
        Returns k random medoids.
        """
        random.sample(range(self.condensed_matrix.row_length),k)
        return random.sample(range(self.condensed_matrix.row_length),k)
        
    def seeding(self,k,seeding_max_cutoff):
        """
        The first k medoids are selected running the gromos algorithm until we get at least
        k clusters.
        """
        print "Seeding process"
        current_cutoff = seeding_max_cutoff
        clusters = []
        while len(clusters)< k and current_cutoff > 0:
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
            self.gromos_clusters_bookkeeping[current_cutoff] = len(clusters)
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
        
