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
 normalize_metrics, score_cluster, get_best_clustering, find_most_negative_cluster

#export PYTHONPATH=/home/victor/git/pyScheduler:$PYTHONPATH

NATIVE_FILE="native.pdb"
TRAJ_PLUS_NATIVE = "traj_plus_native.pdb"
ALL_PLOT = "b_energy_vs_rmsd_all.svg"
REP_50_PLOT = "b_energy_vs_rmsd_50_repr.svg"
BENERGY_RMSD_MATRIX = "be_matrix"

RMSD_COMPRESSION_SCRIPT = "script_compression_dist.json"
CLUSTERING_BE_RMSD_SCRIPT = "script_clustering_be_rmsd.json"
CLUSTERING_TO_5_SCRIPT = "script_cluster_to_5.json"

data = [
#     {
#         'dir':'GPCR-1',
#     },
#     {
#         'dir':'1ewk_1329',
#     },
#     {
#         'dir':'1ewv_1330',
#     },
#     {
#         'dir':'3p0g_1159',
#     },
#      {
#          'dir':'3oe9',
#      },
#     {
#         'dir':'4gpo',
#     },
#     {
#         'dir':'4k5y',
#     },
    {
        'dir':'M1_ACH',
    },
    {
        'dir':'M1_GAM',
    },
    {
        'dir':'M1_NCZ',
    },
    {
        'dir':'M1_XNM',
    }
]

# metrics_to_be_extracted = ["LigandRMSD","L1BindingEne"]
metrics_to_be_extracted = ["L1(63.416.011.9)","L1BindingEne"]


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
    records = []
    processFile(traj_pdb, records, True)
    all_metrics = genMetrics(metrics_to_be_extracted, records)
    print "* %d metrics extracted"%(len(all_metrics))
    plot_metrics(os.path.join(base_dir, ALL_PLOT), all_metrics)

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

    ###########################################################################
    # Do a first rmsd - based compression in order to get the "simplified" plot
    ###########################################################################
    print "* Distance-BASED compression"
    compression_script_path = os.path.join(base_dir, 'scripts', RMSD_COMPRESSION_SCRIPT)
    working_directory = os.path.join(base_dir,"rmsd_compression")
    params = load_dic_in_json(compression_script_path)
    params['global']['pdbs'] = [os.path.join(os.getcwd(), traj_pdb)]
    params['workspace']['base'] = working_directory
    params['matrix']['parameters']['dist_fit_selection'] = "name CA"
    params['matrix']['parameters']['body_selection'] = ligand_description
    save_dic_in_json(params, compression_script_path)
    use_pyproct(working_directory, compression_script_path)

    ############################################################
    # Metrics plotting for the compressed trajectory (50 points)
    ############################################################
    records = []
    processFile(os.path.join(working_directory,"results","compressed_pdb.pdb"), records, True)
    metrics = genMetrics(metrics_to_be_extracted, records)
    print "* %d metrics extracted"%(len(metrics))
    plot_metrics(os.path.join(base_dir, REP_50_PLOT), metrics)

    #######################################################################################################################
    # Generate matrix with metrics (so now we are going to cluster based on Binding Energy and RMSD of the ligand to native
    #######################################################################################################################
    print "* Creating RMSD - Binding Energy matrix"
    matrix_data = scipy.spatial.distance.pdist(normalize_metrics(all_metrics), 'euclidean') # Metrics need to be scaled
    m_handler = MatrixHandler()
    m_handler.distance_matrix = CondensedMatrix(matrix_data)
    matrix_file = os.path.join(base_dir, BENERGY_RMSD_MATRIX)
    m_handler.saveMatrix(matrix_file)

    #######################################################################################################################
    # Cluster by energy-rmsd
    #######################################################################################################################
    print "* RMSD - B. ENERGY clustering"
    be_rmsd_clustering_script_path = os.path.join(base_dir, 'scripts', CLUSTERING_BE_RMSD_SCRIPT)
    last_working_directory = working_directory
    working_directory = os.path.join(base_dir, "benergy_rmsd_clustering")
    params = load_dic_in_json(be_rmsd_clustering_script_path)
    params['global']['pdbs'] = [os.path.join(os.getcwd(), traj_pdb)]
    params['workspace']['base'] = working_directory
    params['matrix']['parameters']['path'] = matrix_file
    save_dic_in_json(params, be_rmsd_clustering_script_path)
    use_pyproct(working_directory, be_rmsd_clustering_script_path)

    #######################################################################################################################
    # Get best clustering
    #######################################################################################################################
    best_clustering = get_best_clustering(working_directory)
    best_clustering = Clustering.from_dic(best_clustering["clustering"])
    scores = []
    for a_cluster in best_clustering.clusters:
        scores.append((score_cluster(a_cluster, m_handler.distance_matrix, all_metrics), a_cluster))

    most_negative = find_most_negative_cluster(scores)
    most_negative_cluster = most_negative[1]

    #######################################################################################################################
    # Plot metrics of most_negative
    #######################################################################################################################
    plot_clusters(os.path.join(base_dir, "clusters.svg"), all_metrics, scores, scores.index(most_negative))

    #######################################################################################################################
    # Store the elements of the most negative cluster
    #######################################################################################################################
    merged_trajectory = os.path.join(working_directory,"tmp","tmp_merged_trajectory.pdb")
    in_handler = open(merged_trajectory,"r")
    most_negative_path = os.path.join(base_dir,"most_negative_cluster.pdb")
    out_handler = open(most_negative_path,"w")
    extract_frames_from_trajectory_sequentially(file_handler_in = in_handler,
                                               number_of_frames = get_number_of_frames(merged_trajectory),
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
    params['global']['pdbs'] = [most_negative_path]
    params['workspace']['base'] = working_directory
    params['matrix']['parameters']['dist_fit_selection'] = "name CA"
    params['matrix']['parameters']['body_selection'] = ligand_description
    save_dic_in_json(params, cluster_to_5_script_path)
    use_pyproct(working_directory, cluster_to_5_script_path)

    #######################################################################################################################
    # Copy the result to the base folder
    #######################################################################################################################
    print "* Copying representatives to base_dir."
    shutil.copyfile(os.path.join(working_directory, "results", "representatives.pdb"), os.path.join(base_dir, "5_representatives.pdb" ))


