"""
Created on 21/05/2012

@author: victor
"""

def find_cutoff_limit(starting_cutoff, min_clusters, max_clusters, grain, hie_algorithm):
    """
    Calculates one of the limits of the cutoff range by adding 'grain' to the starting cutoff until the
    number of clusters of the produced clustering is into the expected range. The direction can be controlled
    with the sign of 'grain'.

    @param starting_cutoff: A guess of the smaller cutoff we can use. It can be 0
    @param min_clusters: The minimum number of clusters allowed by clustering.
    @param max_clusters: The maximum number of clusters allowed by clustering.
    @param grain: The cutoff step that will be used to discover the cutoff range. If it is positive it is going
    to move to the right, if negative to the left. Exiting conditions are also based on this.
    @param hie_algorithm: The hierarchical algorithm instance (holds the hie matrix, so that it has not to be
    calculated each time we get a clustering).
    """
    current_cutoff = starting_cutoff
    search_ended = False
    while not search_ended :
        clustering = hie_algorithm.perform_clustering(kwargs={"cutoff":current_cutoff})
        clustering_size = len(clustering.clusters)
        # Stop when it is into the allowed range
        im_in_the_range = clustering_size >= min_clusters and clustering_size <= max_clusters

        im_out_of_range = False
        if grain>0: # we are moving to the 'right', so to bigger cutoffs and smaller clusters
            im_out_of_range = clustering_size <= min_clusters
        else:
            im_out_of_range = clustering_size >= max_clusters

        search_ended =  im_in_the_range or im_out_of_range
        #print "current", current_cutoff, "->", clustering_size, "   ", grain
        current_cutoff += grain
    return current_cutoff - grain

def get_cutoff_range(starting_cutoff, ending_cutoff, min_clusters, max_clusters, grain, hie_algorithm):
    """
    Returns the cutoff range where the cutoff is inside the range defined by [min_clusters, max_clusters]. It
    first calculates the left limit (from the minimum value of the cutoff, and the smaller number of clusters) and
    then the right limit (from the maximum value of the cutoff).

    @param starting_cutoff: A guess of the smaller cutoff we can use. It can be 0
    @param ending_cutoff: A guess of the bigger cutoff we can use.
    @param min_clusters: The minimum number of clusters allowed by clustering.
    @param max_clusters: The maximum number of clusters allowed by clustering.
    @param grain: Positive number that defines the cutoff step that will be used to discover the cutoff range.
    @param hie_algorithm: The hierarchical algorithm instance (holds the hie matrix, so that it has not to be
    calculated each time we get a clustering).
    """

    lefmost_limit = find_cutoff_limit(starting_cutoff = starting_cutoff,
                                      min_clusters = min_clusters,
                                      max_clusters = max_clusters,
                                      grain = grain,
                                      hie_algorithm = hie_algorithm)

    rightmost_limit = find_cutoff_limit(starting_cutoff = ending_cutoff,
                                          min_clusters = min_clusters,
                                          max_clusters = max_clusters,
                                          grain = -grain,
                                          hie_algorithm = hie_algorithm)

    if lefmost_limit <= rightmost_limit:
        return (lefmost_limit, rightmost_limit)
    else:
        return (rightmost_limit, rightmost_limit+grain)

def get_clusters_with_ranged_search(   hie_algorithm,
                                       cutoff_range_begin,
                                       cutoff_range_end,
                                       min_clusters,
                                       max_clusters,
                                       refine_grain):
    """
    Searches for the range of cutoffs where we are going to get a number of clusters inside the allowed number
    of clusters [min_clusters,max_clusters]. Returns a dictionary indexed by number of clusters where each item
    is a tuple of (cutoff,Clustering). A dictionary is used because all clusterings with same number of clusters
    would be the same (and in this way we get the set without repetition).

    @param cutoff_range_begin: A guess of the smaller cutoff we can use. It can be 0
    @param cutoff_range_end: A guess of the bigger cutoff we can use.
    @param refine_grain: Positive number that defines the cutoff step that will be used to discover the cutoff range.
    @param hie_algorithm: The hierarchical algorithm instance (holds the hie matrix, so that it has not to be
    calculated each time we get a clustering).
    @param min_clusters: The minimum number of clusters allowed by clustering.
    @param max_clusters: The maximum number of clusters allowed by clustering.
    """

    #print cutoff_range_begin, cutoff_range_end,min_clusters,max_clusters
    lefmost_limit, rightmost_limit = get_cutoff_range(cutoff_range_begin,
                                                      cutoff_range_end,
                                                      min_clusters,
                                                      max_clusters,
                                                      0.01,
                                                      hie_algorithm)

    increment = (rightmost_limit-lefmost_limit)/refine_grain

    clusters = {}
    for i in range(refine_grain+1):
        cutoff = increment*i + lefmost_limit
        clustering = hie_algorithm.perform_clustering(kwargs={"cutoff":cutoff})
        clustering_size = len(clustering.clusters)
        if(len(clustering.clusters) != 1):
            clusters[clustering_size]= (cutoff,clustering)

    return clusters