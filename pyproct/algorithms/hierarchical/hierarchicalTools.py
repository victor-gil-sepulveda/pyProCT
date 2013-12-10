'''
Created on 21/05/2012

@author: victor
'''

def get_hierarchical_clustering_for_cutoffs(results,cutoff_list,hie_algorithm,linkage_matrix):
    """
    Does the hierarchical clustering for a list of cutoffs and returns a structure containing the results.
    """
    final_boundary = 0
    for cutoff in cutoff_list:
        clustering = hie_algorithm.perform_clustering(kwargs={"cutoff":cutoff, "hie_mat":linkage_matrix})
        if len(clustering.clusters)>1:
            if final_boundary == 0:
                final_boundary = cutoff
            results.append((cutoff,len(clustering.clusters),clustering))
    return final_boundary


def dicotomic(cutoff_range_begin,
              cutoff_range_end,
              current_cutoff,
              hie_algorithm,
              min_clusters,
              max_clusters):
    """
    Another dicotomic search to get a range of cutoffs where the number of clusters we get from
    the hierarchical clustering is inside the allowed number of clusters range.
    """
    clustering = hie_algorithm.perform_clustering(kwargs={"cutoff":current_cutoff})
    clustering_size = len(clustering.clusters)

    if clustering_size >= min_clusters and clustering_size <= max_clusters:
        return (cutoff_range_begin,cutoff_range_end)

    if clustering_size > max_clusters:
        return dicotomic(current_cutoff,
                         cutoff_range_end,
                         (current_cutoff+cutoff_range_end) / 2,
                         hie_algorithm,
                         min_clusters,
                         max_clusters)

    if clustering_size < min_clusters:
        return dicotomic(cutoff_range_begin,
                         current_cutoff,
                         (current_cutoff+cutoff_range_begin) /2,
                         hie_algorithm,
                         min_clusters,
                         max_clusters)

def get_clusters_with_dicotomic_search(condensed_distance_matrix,
                                       hie_algorithm,
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
    """
    try:
        mrange = dicotomic(cutoff_range_begin,
                           cutoff_range_end,
                           (cutoff_range_begin+cutoff_range_end)/2.,
                           hie_algorithm,
                           min_clusters,
                           max_clusters)
    except RuntimeError:
        print "[WARNING::get_clusters_with_dicotomic_search] dicotomic search needed too much recursion."
        mrange = (cutoff_range_begin,cutoff_range_end)

    _diff = mrange[1]-mrange[0]
    inc = _diff/refine_grain
    clusters = {}
    for i in range(refine_grain+1):
        cutoff = inc*i + mrange[0]
        clustering = hie_algorithm.perform_clustering(kwargs={"cutoff":cutoff})
        clustering_size = len(clustering.clusters)
        if(len(clustering.clusters) != 1):
            clusters[clustering_size]= (cutoff,clustering)
    return clusters