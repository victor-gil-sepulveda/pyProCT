
# http://stackoverflow.com/questions/5147112/matplotlib-how-to-put-individual-tags-for-a-scatter-plot
# http://stackoverflow.com/questions/10374930/matplotlib-annotating-a-3d-scatter-plot

import sys
from pyproclust.clustering.cluster import Cluster
import prody
import matplotlib.pyplot as plt
import numpy
import json
from mpl_toolkits.mplot3d import Axes3D
from pyproclust.tools.commonTools import convert_to_utf8
from pyproclust.clustering.clustering import Clustering
import optparse
import matplotlib.cm as cm

def print_cluster_info(clustering_id, results):
    print clustering_id, results["selected"][clustering_id]["clustering"]["number_of_clusters"], \
    results["selected"][clustering_id]["type"], "".join([ (str(results["scores"][criteria][clustering_id])+", ") for criteria in results["scores"].keys()])

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog -m <arg> -c <arglist> [-o <arg>]', version='1.0')
        
    parser.add_option('-l', action="store_true", dest = "print_list", help="Print a list of generated clusters and some properties",  metavar = "1")
    parser.add_option("-s", action="store_true", dest="show_protein", help="Shows the protein backbone.")
    parser.add_option('-r', action="store", dest = "results_file", help="", metavar = "results.json")
    parser.add_option('-p', action="store", dest = "parameters_file", help="",metavar = "parameters.json")
    parser.add_option('-c', action="store", dest = "clustering_to_see", help="",metavar = "clustering_0001")
    parser.add_option('--stride', type = "int", action="store", dest = "stride", help="",metavar = "5")
    options, args = parser.parse_args()
    
    params = convert_to_utf8(json.loads(open(options.parameters_file).read()))
    results = convert_to_utf8(json.loads(open(options.results_file).read()))
    
    if options.print_list:
        for selected_cluster in results["selected"]:
            print_cluster_info(selected_cluster,results)
        exit()
    
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # Plot protein
    pdb = prody.parsePDB(params["global"]["pdbs"][0])
    if options.show_protein:
        pdb_backbone = pdb.select("name CA").getCoordsets()[0] # "backbone not hetero"
        ax.plot(pdb_backbone.T[0], pdb_backbone.T[1], pdb_backbone.T[2])
    
    # Get geometric centers and plot ligands
    ligand_coords = pdb.select(params["matrix"]["parameters"]["body_selection"]).getCoordsets()
    
    # Get clustering
    if options.clustering_to_see is None:
        options.clustering_to_see = results["best_clustering"]
    clustering = Clustering.from_dic(results["selected"][options.clustering_to_see]["clustering"])
    
    # Print some info
    print_cluster_info(options.clustering_to_see, results)
    
    
    # Show all clusters
    colors = iter(cm.rainbow(numpy.linspace(0, 1, len(clustering.clusters))))
    for cluster in clustering.clusters:
        centers = []
        for i,element in enumerate(cluster.all_elements):
            if options.stride is None or i%options.stride == 0:
                coords = ligand_coords[element]
                centers.append(coords.mean(0))
            
        centers = numpy.array(centers)
        ax.scatter(centers.T[0],centers.T[1],centers.T[2],color=next(colors))
        
    # Plot prototypes
    centers = numpy.array([ligand_coords[cluster.prototype].mean(0) for cluster in clustering.clusters])
    ax.scatter([centers.T[0]],[centers.T[1]],[centers.T[2]], s = 100, c="red", marker='o')
    
    plt.show()
