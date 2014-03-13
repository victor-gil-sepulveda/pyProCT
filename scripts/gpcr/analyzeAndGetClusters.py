import os
import os.path
from scripts.plop.metrics.metrics import processFile, genMetrics
from scripts.gpcr.plotting import plot_metrics, plot_clusters
from pyproct.clustering.clustering import Clustering
from pyproct.tools.pdbTools import extract_frames_from_trajectory_sequentially, get_number_of_frames
from pyRMSD.condensedMatrix import CondensedMatrix
from pyRMSD.matrixHandler import MatrixHandler
import scipy.spatial.distance
import shutil

from scripts.gpcr.tools import load_dic_in_json, save_dic_in_json, use_pyproct,\
 normalize_metrics, score_cluster, get_best_clustering, find_most_negative_cluster,\
    find_5_clusters_with_less_energy


TENERGY_SPAWN_MATRIX = "te_spawn_matrix"
CLUSTERING_SPAWN_TOTE_SCRIPT = "script_clustering_spawn_tenergy.json"
CLUSTERING_TO_5_SCRIPT = "script_cluster_to_5.json"
TOTALE_SPAWN_WORKSPACE = "clustering_totale_spawn"

data = [
#     {
#         'dir':'1ewk',
#         'plots': {
#              "totale_spawning":["L1(32.5-5.933.7)", "TOTALE"],
#              "binding_spawning":["L1(32.5-5.933.7)", "L1BindingEne"]
#          }
#     },
    {
        'dir':'3p0g',
        'plots': {
             "totale_spawning":["L1(63.416.011.9)", "TOTALE"],
             "binding_spawning":["L1(63.416.011.9)", "L1BindingEne"]
         }
    },
    {
        'dir':'3oe9',
        'plots': {
             "totale_spawning":["L1(63.416.011.9)", "TOTALE"],
             "binding_spawning":["L1(63.416.011.9)", "L1BindingEne"]
         }
    },
    {
        'dir':'4gpo',
        'plots': {
             "totale_spawning":["L1(63.416.011.9)", "TOTALE"],
             "binding_spawning":["L1(63.416.011.9)", "L1BindingEne"]
         }
    },
#     {
#         'dir':'4k5y',
#         'plots': {
#              "totale_spawning":["L1(63.416.011.9)", "TOTALE"],
#              "binding_spawning":["L1(63.416.011.9)", "L1BindingEne"]
#          }
#     }
]

cwd = os.getcwd()
for datum in data:
    prot_name = datum['dir']
    print "========================\nWorking with %s\n========================"%(prot_name)
    # Look for the directory and enter it
    base_dir = os.path.join(cwd, prot_name)
    os.chdir(base_dir)
    traj_pdb = "%s.pdb"%prot_name

    ###############################################
    # Generate metrics plot for the whole ensemble
    ###############################################
    plots = datum["plots"]
    for metrics_key in datum["plots"]:
        records = []
        processFile(traj_pdb, records, True)
        all_metrics = genMetrics(plots[metrics_key], records)
        print "* %d metrics extracted"%(len(all_metrics))
        plot_metrics(os.path.join(base_dir, "%s.svg"%metrics_key), all_metrics, plots[metrics_key][0],plots[metrics_key][1])

    ######################################
    # Read the ligand info (only one line)
    ######################################
    ligand_file = open(os.path.join(base_dir, "%s_ligand.txt"%(prot_name)),"r")
    ligand_file_description = ligand_file.readlines()[0].split()
    ligand_file.close()
    ligand = {
              "resname":ligand_file_description[0],
              "chain":ligand_file_description[1],
              "atoms":ligand_file_description[2:]
              }
    ligand_description = "resname %s and name %s"%(ligand["resname"],"".join( a+" " for a in ligand["atoms"]))
    print "* Ligand parsed: ",ligand_description


    #######################################################################################################################
    # Generate matrix with metrics (so now we are going to cluster based on Energy and spawning
    #######################################################################################################################
    print "* Creating Spawning - totalE matrix"
    records = []
    processFile(traj_pdb, records, True)
    all_metrics = genMetrics(plots["totale_spawning"], records)
    matrix_data = scipy.spatial.distance.pdist(normalize_metrics(all_metrics), 'euclidean')
    m_handler = MatrixHandler()
    m_handler.distance_matrix = CondensedMatrix(matrix_data)
    matrix_file = os.path.join(base_dir, TENERGY_SPAWN_MATRIX)
    m_handler.saveMatrix(matrix_file)

    #######################################################################################################################
    # Cluster by metrics
    #######################################################################################################################
    print "* Spawning - totalE clustering"
    be_rmsd_clustering_script_path = os.path.join(base_dir, 'scripts', CLUSTERING_SPAWN_TOTE_SCRIPT)
    working_directory = os.path.join(base_dir, TOTALE_SPAWN_WORKSPACE)
    params = load_dic_in_json(be_rmsd_clustering_script_path)
    params['global']['workspace']['base'] = working_directory
    params['data']['files'] = [os.path.join(os.getcwd(), traj_pdb)]
    params['data']['matrix']['parameters']['path'] = matrix_file
    save_dic_in_json(params, be_rmsd_clustering_script_path)
    use_pyproct(working_directory, be_rmsd_clustering_script_path)

    #######################################################################################################################
    # Get 5 representatives. 2 strategies.
    #######################################################################################################################

    #####################################################################
    #####################################################################
    # Work only with best clustering/ find best cluster
    #####################################################################
    #####################################################################
    best_clustering = get_best_clustering(working_directory)
    best_clustering = Clustering.from_dic(best_clustering["clustering"])
    scores = []
    for a_cluster in best_clustering.clusters:
        scores.append((score_cluster(a_cluster, m_handler.distance_matrix, all_metrics), a_cluster))

    # Remember: first metric is spawning point,  second metric must be energy
    most_negative = find_most_negative_cluster(scores)
    most_negative_cluster = most_negative[1]

    #######################################################################################################################
    # Plot clusters
    #######################################################################################################################
    plot_clusters(os.path.join(base_dir, "clusters.svg"), all_metrics, scores, scores.index(most_negative))

    #######################################################################################################################
    # Store the elements of the most negative cluster
    #######################################################################################################################
    # This works because we have only one traj. Merging before should be mandatory.
    trajectory_path = os.path.join(os.getcwd(), traj_pdb)
    in_handler = open(trajectory_path,"r")
    most_negative_path = os.path.join(base_dir,"most_negative_cluster.pdb")
    out_handler = open(most_negative_path,"w")
    extract_frames_from_trajectory_sequentially(file_handler_in = in_handler,
                                               number_of_frames = get_number_of_frames(trajectory_path),
                                               file_handler_out = out_handler,
                                               frames_to_save = most_negative_cluster.all_elements,
                                               keep_header = True,
                                               write_frame_number_instead_of_correlative_model_number = True)
    in_handler.close()
    out_handler.close()

    #######################################################################################################################
    # Pick 5 well-separated representatives
    #######################################################################################################################
    working_directory = os.path.join(base_dir, "cluster_to_5")
    cluster_to_5_script_path = os.path.join(base_dir, 'scripts', CLUSTERING_TO_5_SCRIPT)
    params = load_dic_in_json(cluster_to_5_script_path)
    params['global']['workspace']['base'] = working_directory
    params['data']['files'] = [most_negative_path]
    params['data']['matrix']['parameters']['dist_fit_selection'] = "name CA"
    params['data']['matrix']['parameters']['body_selection'] = ligand_description
    save_dic_in_json(params, cluster_to_5_script_path)
    use_pyproct(working_directory, cluster_to_5_script_path)

    #######################################################################################################################
    # Copy the result to the base folder
    #######################################################################################################################
    print "* Copying representatives to base_dir."
    shutil.copyfile(os.path.join(working_directory, "results", "representatives.pdb"), os.path.join(base_dir, "5_confs_method_1.pdb" ))

    #####################################################################
    #####################################################################
    # Work with the medoid of the 5 clusterings with lower mean energy
    #####################################################################
    #####################################################################
    print "* Finding medoids."
    clusters = find_5_clusters_with_less_energy(scores)
    medoids = [cluster.calculate_medoid(m_handler.distance_matrix) for en, cluster in clusters]

    #######################################################################################################################
    # Store medoids
    #######################################################################################################################
    print "* Extracting medoids."
    # This works because we have only one traj. Merging before should be mandatory.
    trajectory_path = os.path.join(os.getcwd(), traj_pdb)
    in_handler = open(trajectory_path,"r")
    most_negative_clusters_path = os.path.join(base_dir,"most_negative_clusters.pdb")
    out_handler = open(most_negative_clusters_path,"w")
    extract_frames_from_trajectory_sequentially(file_handler_in = in_handler,
                                               number_of_frames = get_number_of_frames(trajectory_path),
                                               file_handler_out = out_handler,
                                               frames_to_save = medoids,
                                               keep_header = True,
                                               write_frame_number_instead_of_correlative_model_number = True)
    in_handler.close()
    out_handler.close()

    print "* Copying cluster medoids to base_dir."
    shutil.copyfile(most_negative_clusters_path, os.path.join(base_dir, "5_confs_method_2.pdb" ))




