"""
Created on 14/02/2012

@author: victor
"""

"""
Separated static functions of the gromos algorithm for the sake of testing.
"""
def eliminate_cluster_from_node_list(nodes,cluster):
    """
    The representation of the cluster is a tuple with the "center" of the cluster and the list
    of all elements of the cluster.
    """    
    for d in cluster.all_elements:
        nodes.remove(d)