'''
Created on 21/05/2012

@author: victor
'''

def find_cutoff_limit(starting_cutoff, min_clusters, max_clusters, grain, hie_algorithm):
    current_cutoff = starting_cutoff
    search_ended = False
    while not search_ended:
        clustering = hie_algorithm.perform_clustering(kwargs={"cutoff":current_cutoff})
        clustering_size = len(clustering.clusters)
        # Stop when it is into the allowed range
        search_ended =  clustering_size >= min_clusters and clustering_size <= max_clusters
        current_cutoff += grain
    return current_cutoff

def get_cutoff_range(starting_cutoff, ending_cutoff, min_clusters, max_clusters, grain, hie_algorithm):

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

    return (lefmost_limit, rightmost_limit)

def get_clusters_with_ranged_search(   hie_algorithm,
                                       cutoff_range_begin,
                                       cutoff_range_end,
                                       min_clusters,
                                       max_clusters,
                                       refine_grain):
    """
    Searchs for the range of cutoffs where we are going to get a number of clusters inside the allowed number
    of clusters [min_clusters,max_clusters]. Returns a dictionary indexed by number of clusters where each item
    is a tuple of (cutoff,Clustering). A dictionary is used because all clusterings with same number of clusters
    would be the same (and in this way we get the set without repetition).

    @param cutoff_range_begin: A guess of the smaller cutoff we can use. It can be 0
    @param cutoff_range_end: A guess of the bigger cutoff we can use.
    """

    lefmost_limit, rightmost_limit = get_cutoff_range(cutoff_range_begin,
                                                      cutoff_range_end,
                                                      min_clusters,
                                                      max_clusters,
                                                      0.001,
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