'''
Created on 21/05/2012

@author: victor
'''
def calculate_coarse_grain_cutofffs(grain, max_dist):
    """
    Generates a cutoff list dividing the range from 0 to max_dist in grain parts.
    """
    cutoffs = []
    increment =  max_dist / float(grain)
    for i in range(1,grain+1):
        cutoffs.append(max(0,max_dist - i*increment))
    return cutoffs
  
def calculate_fine_grain_cutoffs(coarse_grain,fine_grain,max_dist, final_boundary):
    """
    Does the same division but for a subrange.
    """
    search_range = []
    increment =  max_dist / float(coarse_grain)
    fine_increment = increment/fine_grain
    for i in range(1,fine_grain):
        search_range.append(final_boundary+(fine_increment*i))
    search_range.reverse()
    return search_range

def dico(search_range,i,j,k,hie_algorithm,linkage_matrix):
    """
    Dicotomic search for an ordered list of cutoffs. It shows a more exact point where the
    number of clusters starts to be > 1
    """
    clustering = hie_algorithm.perform_clustering(kwargs={"cutoff":search_range[k], "hie_mat":linkage_matrix})
    clustering_size = len(clustering.clusters)
    if i<j-2:
        if clustering_size == 1: # Then we have to look in the [k, j] (right) zone because of the ordering
            return dico(search_range,k,j,int((j+k)/2.),hie_algorithm,linkage_matrix)
        else: # We search in the left half ( [i,k] zone)
            # We keep picking the results
            return dico(search_range,i,k,int((i+k)/2.),hie_algorithm,linkage_matrix)
    else:
        return k

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

def calculate_cutoff_numcluster_list(hie_algorithm, linkage_matrix, max_dist, coarse_grain = 100, fine_grain = 100):
    """
    Calculates the clusterings for different lists of cutoffs, searching for the point
    where the number of clusters becomes larger than 1.
    Once this point has been found, it makes a deeper search inside the range defined by
    the previous cutoff and the next.
    For instance if we have a max dist of 6 and coarse grain of 4, the cutoffs will be:
    [0, 1.5, 3, 4.5]
    If the point where we have more than 2 cluster is 1.5, we'll perform a fine grain search
    between 1.5 and 3.
    """
    results = []
    # Coarse grain search
    coarse_grain_cutoffs = calculate_coarse_grain_cutofffs(coarse_grain, max_dist)
    final_boundary = get_hierarchical_clustering_for_cutoffs(results,coarse_grain_cutoffs, hie_algorithm, linkage_matrix)
    # Now refine the search
    fine_grain_cutoffs =  calculate_fine_grain_cutoffs(coarse_grain,fine_grain,max_dist, final_boundary)
    # Better doing a dicotomic-like search
    bound_start = dico(fine_grain_cutoffs,0,fine_grain,int(fine_grain/2),hie_algorithm,linkage_matrix)
    interesting_cutoffs = fine_grain_cutoffs[bound_start:]
    get_hierarchical_clustering_for_cutoffs(results,interesting_cutoffs, hie_algorithm, linkage_matrix)
    return results


def dicotomic(cutoff_range_begin,cutoff_range_end, current_cutoff, hie_algorithm,min_clusters,max_clusters):
    """
    Another dicotomic search to get a range of cutoffs where the number of clusters we get from
    the hierarchical clustering is inside the allowed number of clusters range.
    """
    clustering = hie_algorithm.perform_clustering(kwargs={"cutoff":current_cutoff})
    clustering_size = len(clustering.clusters)
    
    if clustering_size >= min_clusters and clustering_size <= max_clusters:
        return (cutoff_range_begin,cutoff_range_end)
    
    if clustering_size > max_clusters:
        return dicotomic(current_cutoff,cutoff_range_end, (current_cutoff+cutoff_range_end) / 2, hie_algorithm,min_clusters,max_clusters)
    
    if clustering_size < min_clusters:
        return dicotomic(cutoff_range_begin,current_cutoff, (current_cutoff+cutoff_range_begin) /2, hie_algorithm,min_clusters,max_clusters)
        
def get_clusters_with_dicotomic_search(condensed_distance_matrix,hie_algorithm,cutoff_range_begin,cutoff_range_end,min_clusters,max_clusters,refine_grain):
    """
    Searchs for the range of cutoffs where we are going to get a number of clusters inside the allowed number
    of clusters [min_clusters,max_clusters]. Returns a dictionary indexed by number of clusters where each item 
    is a tuple of (cutoff,Clustering). A dictionary is used because all clusterings with same number of clusters 
    would be the same (and in this way we get the set without repetition). 
    """
    mrange = dicotomic(cutoff_range_begin, cutoff_range_end, (cutoff_range_begin+cutoff_range_end)/2., hie_algorithm, min_clusters, max_clusters)
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