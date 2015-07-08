import numpy
import sys
from pyproct.postprocess.actions.confSpaceComparison.tools import getAllElements
from pyproct.clustering.evaluation.metrics.common import get_distances_of_elements_to

class OverlapCalculator():

    def __init__(self):
        """
        Constructor. This class works as namespace.
        """
        pass

    @classmethod
    def calculate_clustering_overlap(cls, decomposed_clusters, distance_matrix):
        """

        """
        cluster_overlaps = numpy.array([ cls.calculate_cluster_overlap(decomposed_cluster, distance_matrix) for decomposed_cluster in decomposed_clusters])
        cluster_sizes = numpy.array([len(getAllElements(decomposed_cluster)) for decomposed_cluster in decomposed_clusters])

        num_elements = numpy.sum(cluster_sizes)
        clustering_overlap =  ((1./num_elements)* numpy.dot(cluster_overlaps, cluster_sizes))
        
        # Overlap 0 is the best overlap, overlap 1 is the worst. We reverse it to do it more easy to understand
        return 1 - clustering_overlap

    @classmethod
    def calculate_cluster_overlap(cls, decomposed_cluster, distance_matrix):
        """
        Calculates the overlap value for a cluster in a range [0,1].

        @param decomposed_cluster: A
        """
        if len(decomposed_cluster) == 1:
                return 1.0 # If the cluster is 'pure' we penalize the global overlap

        else:
            N = len(getAllElements(decomposed_cluster))
                    
            min_distances = cls.get_cluster_min_distances(decomposed_cluster, distance_matrix)
            
            max_min_distance = max(min_distances)
            
            if max_min_distance == 0: # Then overlap is total
                return 0.0
            else:
                return (1./N) * (numpy.sum(min_distances) /max_min_distance)

    @classmethod
    def get_cluster_min_distances(cls, decomposed_cluster, distance_matrix):
        """
        Calculates the distances between the elements of all different classes in the cluster
        and returns the minimum distance for each of these elements.
        Some distances will be counted twice. This is OK.
        """
        allIds = decomposed_cluster.keys()
        min_distances = []

        if len(allIds)>1: # if the cluster is pure, we do not calculate min or max (it does not have sense)
            for setId in allIds:
                myVsIds = list(allIds)
                myVsIds.remove(setId)
                vs_elements = []
                for vsId in myVsIds:
                    vs_elements.extend(decomposed_cluster[vsId])
                for element in decomposed_cluster[setId]:
                    min_distances.append( numpy.min(get_distances_of_elements_to(element, vs_elements, distance_matrix)))
            return numpy.array(min_distances)
        else:
            raise ValueError("Asking min max distances of a PURE cluster.")

