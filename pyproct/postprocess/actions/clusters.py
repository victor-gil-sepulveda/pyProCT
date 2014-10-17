"""
Created on 16/07/2014

@author: victor
"""

import os.path
import tarfile
from pyproct.postprocess.actions.representatives import save_cluster_elements

class SaveAllClustersPostAction(object):
    KEYWORD = "representatives"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, data_handler, workspaceHandler, matrixHandler, generatedFiles):

        results_place = workspaceHandler["results"]
        clusters_place = workspaceHandler["clusters"]
    
        cluster_files = []
        for cluster in clustering.clusters:
            output_path = os.path.join(clusters_place, 
                                       "%s.pdb"%(cluster.id))
            save_cluster_elements(cluster.all_elements,
                                  [cluster.id]*len(cluster.all_elements), # all share same cluster id
                                  output_path,
                                  data_handler,
                                  postprocessing_parameters)
            cluster_files.append(output_path)
    
        # Add all bz2 files to a tar file
        tar_path = os.path.join(results_place,
                                "clusters.tar.gz")
        
        tar = tarfile.open(tar_path, "w:gz")
        for comp_file in cluster_files:
            tar.add(comp_file, os.path.basename(comp_file))
        tar.close()
    
        generatedFiles.append({
                               "description":"Clusters",
                               "path":os.path.abspath(tar_path),
                               "type":"compressed_pdb"
        })
