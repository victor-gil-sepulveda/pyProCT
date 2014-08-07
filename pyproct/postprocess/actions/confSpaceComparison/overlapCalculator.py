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
    def calculate_global_overlap(cls, decomposed_clusters, distance_matrix, cluster_method, global_method):
        """

        """
        sys.stdout.flush()
        cluster_overlaps = numpy.array([ cls.calculate_cluster_overlap(cluster_method, decomposed_cluster, distance_matrix) for decomposed_cluster in decomposed_clusters])
        cluster_sizes = numpy.array([len(getAllElements(decomposed_cluster)) for decomposed_cluster in decomposed_clusters])
        global_size = numpy.sum(cluster_sizes)
        overlap = numpy.sum(cluster_overlaps * cluster_sizes) / float(global_size)
        # Overlap 0 is the best overlap, overlap 1 is the worst. We reverse it to do it more easy to understand
        return 1 - overlap

    @classmethod
    def calculate_cluster_overlap(cls, method, decomposed_cluster, distance_matrix):
        """
        Calculates the overlap value for a cluster in a range [0,1].

        @param decomposed_cluster: A
        """
        min_distances, max_distances= cls.get_cluster_min_max_distances(decomposed_cluster, distance_matrix)

        if len(min_distances) == 0:
            return 1.

        if method == 1:
            return numpy.sum(min_distances) / numpy.sum(max_distances)

        elif method == 2:
            return numpy.sum(min_distances / max_distances) / len(getAllElements(decomposed_cluster))

        else:
            print "[ERROR OverlapCalculator::calculate_cluster_overlap] The method nr. %d does not exist."%(method)
            exit()

    @classmethod
    def get_cluster_min_max_distances(cls, decomposed_cluster, distance_matrix):
        """

        """
        allIds = decomposed_cluster.keys()
        min_distances = []
        max_distances = []

        if len(allIds)>1: # if the cluster is pure, we do not calculate min or max (it does not have sense)
            for setId in allIds:
                myVsIds = list(allIds)
                myVsIds.remove(setId)
                vs_elements = []
                for vsId in myVsIds:
                    vs_elements.extend(decomposed_cluster[vsId])
                for element in decomposed_cluster[setId]:
                    min_distances.append( numpy.min(get_distances_of_elements_to(element, vs_elements, distance_matrix)))
                    max_distances.append( numpy.max(get_distances_of_elements_to(element, vs_elements, distance_matrix)))
        return numpy.array(min_distances), numpy.array(max_distances)



