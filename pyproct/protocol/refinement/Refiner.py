'''
Created on 11/02/2013

@author: victor
'''

#             ######################
#             # Refine results
#             ######################
#             if protocol_params.do_refinement:
#                 ## Change params
#                 old_min = protocol_params.min_clusters
#                 old_max = protocol_params.max_clusters
#                 old_kmedoids_step = protocol_params.kmedoids_step
#                 old_spectral_step = protocol_params.spectral_clustering_step
#                 protocol_params.max_clusters = protocol_params.refinement_max_clusters
#                 protocol_params.kmedoids_step = protocol_params.refinement_step
#                 protocol_params.spectral_clustering_step = protocol_params.refinement_step
#                 
#                 separator = Separator()
#                 pure_A , pure_B, mixed_with_elements = separator.separate(best_clustering, protocol_params.pdb1, protocol_params.pdb2)
#                 pre_refinement_cluster_lengths = {}
#                 pre_refinement_cluster_lengths['pure_A'] = len(pure_A)
#                 pre_refinement_cluster_lengths['pure_B'] = len(pure_B)
#                 pre_refinement_cluster_lengths['mixed'] = len(mixed_with_elements)
#                 
#                 print "[refinement] Initial params, max_clusters: ", protocol_params.max_clusters,"step: ",protocol_params.kmedoids_step
#                 print "[refinement] Initial clusters A: %d B: %d Mixed:%d"%(len(pure_A),len(pure_B),len(mixed_with_elements)) 
#                 
#                 new_clusters = []
#                 post_refinement_cluster_lenghts = {}
#                 
#                 new_pure_A = None
#                 if len(pure_A) != 0:
#                     print "[refinement] Refining A"
#                     protocol_params.min_clusters = max(protocol_params.refinement_min_clusters,len(pure_A))
#                     print "[refinement] Refining A. protocol_params.min_clusters changed from %d to %d"%(old_min,protocol_params.min_clusters)
#                     refiner = pureRefinementProtocol(pure_A, self.matrixHandler.distance_matrix,\
#                                                      self.workspaceHandler.refinement_pure_A+'/clusterings',\
#                                                      self.trajectoryHandler.pdb_structure)
#                     new_pure_A = refiner.run(protocol_params)
#                     if new_pure_A != None:
#                         new_clusters.extend(new_pure_A)
#                         post_refinement_cluster_lenghts['pure_A'] = len(new_pure_A)
#                         print "[refinement] A has been refined. Clusters changed from %d to %d"%(len(pure_A),len(new_pure_A))
#                     else:
#                         print "[refinement] Impossible to refine A"
#                 else:
#                     print "[refinement] No clusters in A, impossible to refine."
#                 
#                 if new_pure_A == None:
#                     new_clusters.extend(pure_A)
#                     post_refinement_cluster_lenghts['pure_A'] = pre_refinement_cluster_lengths['pure_A']
#                 
#                 new_pure_B = None   
#                 if len(pure_B) != 0:
#                     print "[refinement] Refining B"
#                     protocol_params.min_clusters = max(protocol_params.refinement_min_clusters,len(pure_B))
#                     print "[refinement] Refining B. protocol_params.min_clusters changed from %d to %d"%(old_min,protocol_params.min_clusters)
#                     refiner = pureRefinementProtocol(pure_B, self.matrixHandler.distance_matrix,\
#                                                      self.workspaceHandler.refinement_pure_B+'/clusterings',\
#                                                      self.trajectoryHandler.pdb_structure)
#                     new_pure_B = refiner.run(protocol_params)
#                     if new_pure_B != None:
#                         new_clusters.extend(new_pure_B)
#                         post_refinement_cluster_lenghts['pure_B'] = len(pure_B)
#                         print "[refinement] B has been refined. Clusters changed from %d to %d"%(len(pure_B),len(new_pure_B))
#                     else:
#                         print "[refinement] Impossible to refine B"
#                 else:
#                     print "[refinement] No clusters in B, impossible to refine."
#                     
#                 if new_pure_B == None:
#                     new_clusters.extend(pure_B)
#                     post_refinement_cluster_lenghts['pure_B'] = pre_refinement_cluster_lengths['pure_B']
#                 
#                 new_mixed_clusters = None
#                 if len(mixed_with_elements) != 0:
#                     print "[refinement] Refining Mixed"
#                     protocol_params.min_clusters = max(protocol_params.refinement_min_clusters,len(mixed_with_elements))
#                     print "[refinement] Refining Mixed. protocol_params.min_clusters changed from %d to %d"%(old_min,protocol_params.min_clusters)
#                     
#                     refiner = mixedRefinementProtocol(mixed_with_elements, self.matrixHandler.distance_matrix,\
#                                                      self.workspaceHandler.refinement_mixed+'/clusterings',\
#                                                      self.trajectoryHandler.pdb_structure)
#                     new_mixed_clusters = refiner.run(protocol_params)
#                     if new_mixed_clusters != None:
#                         new_pure_A, new_pure_B, new_mixed = new_mixed_clusters
#                         new_clusters.extend(new_mixed)
#                         new_clusters.extend(new_pure_A)
#                         new_clusters.extend(new_pure_B)
#                         post_refinement_cluster_lenghts['mixed'] = len(new_mixed)
#                         post_refinement_cluster_lenghts['mixed_pure_A'] = len(new_pure_A)
#                         post_refinement_cluster_lenghts['mixed_pure_B'] = len(new_pure_B)
#                         print "[refinement] Mixed has been refined. Clusters changed from %d to %d"%(len(mixed_with_elements),len(new_mixed))
#                     else:
#                         print "[refinement] Impossible to refine B"
#                 else:
#                     print "[refinement] No clusters in B, impossible to refine." 
#                   
#                 if new_mixed_clusters == None:
#                     only_mixed = []
#                     for pack in mixed_with_elements:
#                         only_mixed.append(pack[0])
#                     new_clusters.extend(only_mixed )
#                     post_refinement_cluster_lenghts['mixed'] = pre_refinement_cluster_lengths['mixed']
#                     post_refinement_cluster_lenghts['mixed_pure_A'] = 0
#                     post_refinement_cluster_lenghts['mixed_pure_B'] = 0
#                 
#                 open(self.workspaceHandler.results_path+"/refinement.txt","w").write("""pure A: %d -> %d
# pure B: %d -> %d
# mixed: %d -> %d
# new pure A in mixed: %d
# new pure B in mixed: %d\n"""%(pre_refinement_cluster_lengths["pure_A"],post_refinement_cluster_lenghts["pure_A"],\
#                          pre_refinement_cluster_lengths["pure_B"],post_refinement_cluster_lenghts["pure_B"],\
#                          pre_refinement_cluster_lengths["mixed"],post_refinement_cluster_lenghts["mixed"],\
#                          post_refinement_cluster_lenghts['mixed_pure_A'],post_refinement_cluster_lenghts['mixed_pure_B']))
#                 
#                 best_clustering = Clustering(new_clusters, "Refined")
#                 all_elems = []
#                 for c in best_clustering.clusters:
#                     all_elems.extend(c.all_elements)
#                 
#                 ## Restore params
#                 protocol_params.max_clusters = old_max
#                 protocol_params.kmedoids_step = old_kmedoids_step
#                 protocol_params.spectral_clustering_step = old_spectral_step