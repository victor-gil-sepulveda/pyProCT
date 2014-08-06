"""
Created on 16/05/2012

@author: victor
"""
from pygraph.classes.digraph import digraph #@UnresolvedImport
from pygraph.readwrite.dot import write #@UnresolvedImport
import pyproct.clustering.cluster as clusTools
import pyproct.tools.commonTools as common
import os
from pyproct.clustering.cluster import Cluster
from pyproct.tools.distributionTools import get_distance_std_dev_for_elems
from pyproct.clustering.clustering import Clustering

def gen_color(num_elems,max_elems):
    """
    Generates a red color in 'dot' format, which tone is based in the number of elements of a cluster (more elements, more intense).
    
    @param num_elems: 
    @param max_elems: 
    
    @return: 
    """
    red_tone = num_elems / float(max_elems)
    color = "0.000 %.3f 1.000"%(red_tone)
#    print color
    return color.lower()

def populate_nodes_with_labels(clustering, num_elems_of_traj_2, std_deviations, graph):
    """
    """
#    print "num_elems_of_traj_2", num_elems_of_traj_2
    labels = []
    num_clusters = len(clustering.clusters)
    max_elems = max(num_elems_of_traj_2)
    for i in range(num_clusters):
#        print clustering.clusters[i].get_size(), " ",  num_elems_of_traj_2[i]
        label = "cluster_"+str(i)+" ("+str(clustering.clusters[i].get_size())+"/"+str(num_elems_of_traj_2[i])+")"
        labels.append(label)
        size = std_deviations[i]*10+1
        color = gen_color(num_elems_of_traj_2[i],max_elems)
        #print color
        graph.add_node(label,[("shape","circle"),( "style","filled"),("fillcolor",color),("color","black"),("width",str(size))])
    return labels

def calculate_probability_matrix(clustering):
    """
    """
    num_clusters = len(clustering.clusters)
    total_elements,cluster_sizes = clusTools.get_cluster_sizes(clustering.clusters) #@UnusedVariable
    class_list = clustering.gen_class_list()
    
    prob_matrix = []
    for i in range(num_clusters):
        row = [0.]*num_clusters
        prob_matrix.append(row)
    
    prob_increments = []
    for i in range(num_clusters):
        prob_increments.append(1./cluster_sizes[i])
    
    for i in range(len(class_list)-1):
        current_cluster = class_list[i]
        next_cluster = class_list[i+1]
        prob_matrix[current_cluster][next_cluster] += prob_increments[current_cluster]
    
    return prob_matrix

def add_graph_edges(graph,labels,clustering,prob_matrix):
    """
    """
    num_clusters = len(clustering.clusters)
    for i in range(num_clusters):
        for j in range(num_clusters):
            prob = prob_matrix[i][j]
            if prob!= 0 and prob > 0.05:
                graph.add_edge((labels[i],labels[j]),wt=prob,label = '%.2f%%'%(prob*100))

def do_graph(clustering,num_elems_of_traj_2,std_deviations,filename):
    """
    """
    graph = digraph()
    labels = populate_nodes_with_labels(clustering, num_elems_of_traj_2,std_deviations, graph)
    prob_matrix = calculate_probability_matrix(clustering)
    add_graph_edges(graph,labels,clustering,prob_matrix)
    tmp_file = open("tmp_dot","w")
    tmp_file.write(write(graph))
    tmp_file.close()
    common.print_and_flush("delegating to dot...")
    os.system("cat tmp_dot|dot -Tpng -o "+filename+";rm tmp_dot")
 
def purge_mixed_clusters_and_do_graph(mixed, pure_clusters_traj1,condensed_distance_matrix,std_devs_from_A,path):
    """
    """
    common.print_and_flush( "Purging clusters...")
    # Purge all mixed clusters of elements from traj2
    purged = []
    num_elems_of_traj_2 = []
    for i in range(len(mixed)):
        cluster, elems_in_traj1, elems_in_traj2 = mixed[i] #@UnusedVariable
        num_elems_of_traj_2.append(len(elems_in_traj2))
        # We rebuild the cluster with only elements of traj 1
        purged.append(Cluster(prototype=None,elements = elems_in_traj1))
#        print "l ",len(elems_in_traj1)," ",len(elems_in_traj2)
    
    # we also need to have traj 1 pure clusters
    purged.extend(pure_clusters_traj1)
    
    # Those don't have any element of traj 2, so we put 0s in the number of 
    # elements list
    num_elems_of_traj_2.extend([0]*len(pure_clusters_traj1))
    
    #Calculate statistics for the remaining clusters
    for i in range(len(pure_clusters_traj1)):
        medoid = pure_clusters_traj1[i].calculate_medoid(condensed_distance_matrix)
        std_devs_from_A.append(get_distance_std_dev_for_elems(pure_clusters_traj1[i].all_elements,medoid,condensed_distance_matrix))
    common.print_and_flush( "Done.\n")
    
    common.print_and_flush("Trying to draw state graph...")
    do_graph(Clustering(purged,sort =  False),num_elems_of_traj_2,std_devs_from_A,path)
    common.print_and_flush("Done.\n")
