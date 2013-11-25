'''
Created on Mar 26, 2013

@author: victor
'''
import os
from pyproclust.driver.handlers.timerHandler import TimerHandler
from pyproclust.driver.handlers.workspaceHandler import WorkspaceHandler
from pyproclust.driver.handlers.trajectoryHandler import TrajectoryHandler
from pyproclust.driver.handlers.matrix.matrixHandler import MatrixHandler
from pyproclust.driver.observer.observable import Observable
import pyproclust.tools.plotTools as plotTools
from pyproclust.protocol.protocol import ClusteringProtocol
import pyproclust.protocol.saveTools as saveTools
from pyproclust.clustering.comparison.distrprob.kullbackLieblerDivergence import KullbackLeiblerDivergence
from pyproclust.driver.compressor.compressor import Compressor
from pyproclust.driver.results.clusteringResultsGatherer import ClusteringResultsGatherer
import json
from pyproclust.clustering.clustering import Clustering
from pyproclust.clustering.comparison.caDisplacement import CA_mean_square_displacement_of_cluster
from pyproclust.clustering.cluster import Cluster
from pyRMSD.RMSDCalculator import RMSDCalculator
import numpy
import matplotlib.cm as cm

class Driver(Observable):

    def __init__(self, observer):
        super(Driver, self).__init__(observer)
        self.generatedFiles = []

    def create_workspace(self, parameters):
        self.workspaceHandler = WorkspaceHandler(parameters["workspace"], self.observer)
        self.workspaceHandler.create_directories()

    def save_parameters_file(self, parameters):
        parameters_file_path = os.path.join(self.workspaceHandler["results"], "parameters.json")
        open(parameters_file_path, "w").write(json.dumps(parameters.params_dic, sort_keys=False, indent=4, separators=(',', ': ')))
        self.generatedFiles = [{"description":"Parameters file", "path":parameters_file_path, "type":"text"}]

    def create_matrix(self, parameters):
        self.matrixHandler = MatrixHandler(parameters["matrix"])
        self.notify("Matrix calculation", [])
        self.timer.start("Matrix Generation")
        self.matrixHandler.create_matrix(self.trajectoryHandler)
        statistics_file_path = self.matrixHandler.save_statistics(self.workspaceHandler["matrix"])
        self.generatedFiles.append({"description":"Matrix statistics", "path":statistics_file_path, "type":"text"})
        self.timer.stop("Matrix Generation")
        if "filename" in parameters["matrix"]:
            self.timer.start("Matrix Save")
            self.matrixHandler.save_matrix(os.path.join(self.workspaceHandler["matrix"], parameters["matrix"]["filename"]))
            self.timer.stop("Matrix Save")
        #########################
        # Matrix plot
        #########################
        if "image" in parameters["matrix"].keys():
            self.timer.start("Matrix Imaging")
            matrix_image_file_path = os.path.join(self.workspaceHandler["matrix"], parameters["matrix"]["image"]["filename"])
            plotTools.matrixToImage(self.matrixHandler.distance_matrix, matrix_image_file_path, max_dim=parameters["matrix"]["image"]["dimension"], observer=self.observer)
            self.generatedFiles.append({"description":"Matrix image", "path":matrix_image_file_path, "type":"image"})
            self.timer.stop("Matrix Imaging")

    def get_best_clustering(self, parameters):
        best_clustering = None
        ##############################
        # Do the actual clustering
        ##############################
        clustering_results = None

        if parameters["clustering"]["generation"]["method"] == "generate":
            clustering_results = ClusteringProtocol(self.timer, self.observer).run(parameters, self.matrixHandler,
                                                                                                self.workspaceHandler,
                                                                                                self.trajectoryHandler)
            if clustering_results != None:
                best_clustering = None
                best_clustering_id, selected, not_selected, scores = clustering_results
                best_clustering = selected[best_clustering_id]
            else:
                self.notify("SHUTDOWN", "Improductive clustering search. Relax evaluation constraints.")
                print "[FATAL Driver:get_best_clustering] Improductive clustering search. Exiting..."
                exit()

        ##############################
        # Load the clustering
        ##############################
        if parameters["clustering"]["generation"]["method"] == "load":
            best_clustering = {"clustering":Clustering.from_dic(parameters["clustering"]["generation"])}

        #################################
        # Abort if no clusters were found
        #################################
        if best_clustering is None:
            self.notify("SHUTDOWN", "Improductive clustering search. Relax evaluation constraints.")
            print "[FATAL Driver:get_best_clustering] Improductive clustering search. Exiting..."
            exit()

        return best_clustering, clustering_results

    def save_clustering_results(self, clustering_results):
        results_path = os.path.join(self.workspaceHandler["results"], "results.json")
        self.generatedFiles.append({"description":"Results file", "path":results_path, "type":"text"})
        json_results = ClusteringResultsGatherer().gather(self.timer,
                                                            self.trajectoryHandler,
                                                            self.workspaceHandler,
                                                            clustering_results,
                                                            self.generatedFiles)
        # Results are first added and saved later to avoid metareferences :D
        open(results_path, "w").write(json_results)



    def postprocess(self, parameters, best_clustering):
        ##############################
        # Specialized post-processing
        ##############################
        action_type = parameters["global"]["action"]["type"]

        if action_type == "clustering" or action_type == "advanced":
            ##############################
            # Saving representatives
            ##############################
            self.timer.start("Representatives")
            medoids = best_clustering["clustering"].get_medoids(self.matrixHandler.distance_matrix)
            # Set prototypes and ids (medoids are ordered)
            for i in range(len(best_clustering["clustering"].clusters)):
                best_clustering["clustering"].clusters[i].prototype = medoids[i]

            keep_remarks = False
            try:
                keep_remarks = parameters["global"]["action"]["parameters"]["keep_remarks"]
            except:
                pass

            #Saver CA mean squared displacement of best cluster
            #TODO: REFACTORING
#             try:
            if parameters["matrix"]["method"] == "rmsd":
                global_cluster = Cluster(None, best_clustering["clustering"].get_all_clustered_elements())
                global_cluster.prototype = global_cluster.calculate_medoid(self.matrixHandler.distance_matrix)
                ca_pdb_coordsets =numpy.copy(self.trajectoryHandler.getJoinedPDB().select("name CA").getCoordsets())
                calculator = RMSDCalculator(calculatorType = "QTRFIT_SERIAL_CALCULATOR",
                                                fittingCoordsets = ca_pdb_coordsets)
                calculator.iterativeSuperposition()
                CA_mean_square_displacements= {
                                               "global":list(CA_mean_square_displacement_of_cluster(ca_pdb_coordsets,\
                                                                                                    global_cluster))
                                               }
                clusters = best_clustering["clustering"].clusters
                for i in range(len(clusters)):
                    cluster = clusters[i]
                    # Pick the coordinates (ensuring that we are copying them)
                    fitting_coordinates_of_this_cluster = ca_pdb_coordsets[cluster.all_elements]
                    calculator = RMSDCalculator(calculatorType = "QTRFIT_SERIAL_CALCULATOR",
                                                fittingCoordsets = fitting_coordinates_of_this_cluster)

                    # Make an iterative superposition (to get the minimum RMSD of all with respect to a mean conformation)
                    calculator.iterativeSuperposition()

                    # Calculate and convert to list (to serialize)
                    CA_mean_square_displacements[cluster.id] = list(CA_mean_square_displacement_of_cluster(ca_pdb_coordsets,\
                                                                                                           cluster))

                displacements_path = os.path.join(self.workspaceHandler["results"], "CA_displacements.json")

                self.generatedFiles.append({
                                            "description":"Alpha Carbon mean square displacements",
                                            "path":displacements_path,
                                            "type":"text"
                })

                open(displacements_path,"w").write(json.dumps(CA_mean_square_displacements,
                                                      sort_keys=False,
                                                      indent=4,
                                                      separators=(',', ': ')))
#             except Exception:
#                 print "Impossible to calculate CA displacements"

            if parameters["matrix"]["method"] == "distance":
            	# TODO: Superpose again

            	centers_path = os.path.join(self.workspaceHandler["results"], "selection_centers.json")

                self.generatedFiles.append({
                                            "description":"Centers of the selection used to calculate distances",
                                            "path":centers_path,
                                            "type":"text"
                })

                clustering = best_clustering["clustering"]
                ligand_coords = self.trajectoryHandler.getSelection(parameters["matrix"]["parameters"]["body_selection"])

                centers_contents={}
                centers = []
                # Superpose and center coords

                # Center coords
                #for i in range(len(ligand_coords)):
                #    ligand_coords[i] -= ligand_coords[i].mean(0)

                # Get Bounding Box
                [max_x,max_y,max_z] = numpy.max(numpy.max(ligand_coords,1),0)
                [min_x,min_y,min_z] = numpy.min(numpy.min(ligand_coords,1),0)
                centers_contents["bounding_box"] = [  [max_x, max_y, max_z],
                                            [max_x, max_y, min_z],
                                            [max_x, min_y, max_z],
                                            [max_x, min_y, min_z],
                                            [min_x, max_y, max_z],
                                            [min_x, max_y, min_z],
                                            [min_x, min_y, max_z],
                                            [min_x, min_y, min_z]]
#for p1 in ["max_x","min_x"]:
#    for p2 in ["max_y","min_y"]:
#        for p3 in ["max_z","min_z"]:
#            print "[%s, %s, %s],"%(p1,p2,p3)
#
                colors = iter(cm.rainbow(numpy.linspace(0, 1, len(clustering.clusters))))
                centers_contents["points"]
                for cluster in clustering.clusters:
                    centers = []
                    for i,element in enumerate(cluster.all_elements):
                        coords = ligand_coords[element]
                        centers.append(list(coords.mean(0)))
                    centers_contents["points"][cluster.id] = {}
                    centers_contents["points"][cluster.id]["prototype"] = list(ligand_coords[cluster.prototype].mean(0))
                    centers_contents["points"][cluster.id]["centers"] = centers
                    centers_contents["points"][cluster.id]["color"] = list(next(colors))[0:3]

                open(centers_path,"w").write(json.dumps(centers_contents,
                                                      sort_keys=False,
                                                      indent=4,
                                                      separators=(',', ': ')))


            representatives_path = saveTools.save_representatives(medoids,
                                                                  "representatives",
                                                                  self.workspaceHandler,
                                                                  self.trajectoryHandler,
                                                                  do_merged_files_have_correlative_models=True,
                                                                  write_frame_number_instead_of_correlative_model_number=True,
                                                                  keep_remarks = keep_remarks )
            self.generatedFiles.append({"description":"Cluster central conformations",
                                        "path":representatives_path,
                                        "type":"pdb"})
            self.timer.stop("Representatives")

        elif action_type == "comparison":
            ############################################
            # Distribution analysis
            ############################################
            self.timer.start("KL divergence")
            klDiv = KullbackLeiblerDivergence(self.trajectoryHandler.pdbs, self.matrixHandler.distance_matrix)
            kl_file_path = os.path.join(self.workspaceHandler["matrix"], "kullback_liebler_divergence")
            klDiv.save(kl_file_path)
            matrix_image_file_path = os.path.join(self.workspaceHandler["matrix"], parameters["matrix"]["image"]["filename"])
            self.generatedFiles.append({"description":"Kullback-Leibler divergence",
                                        "path":matrix_image_file_path,
                                        "type":"text"})
            self.timer.stop("KL divergence")

        elif action_type == "compression":
            ############################################
            # Compress
            ############################################
            self.timer.start("Compression")
            compressor = Compressor(parameters["global"]["action"]["parameters"])
            compressed_file_path = compressor.compress(best_clustering["clustering"],
                                                       (lambda params: params['file'] if 'file' in params else "compressed_pdb")(parameters["global"]["action"]["parameters"]),
                                                       self.workspaceHandler,
                                                       self.trajectoryHandler,
                                                       self.matrixHandler)
            self.generatedFiles.append({"description":"Compressed file",
                                        "path":compressed_file_path,
                                        "type":"pdb"})
            self.timer.stop("Compression")

    def perform_actions(self, parameters):
        best_clustering, clustering_results = self.get_best_clustering(parameters)
        self.postprocess(parameters, best_clustering)
        #################################
        # Results are saved to a file
        #################################
        self.save_clustering_results(clustering_results)


    def run(self, parameters):

        #####################
        # Start timing
        #####################
        self.timer = TimerHandler()
        self.timer.start("Global")

        #####################
        # Workspace Creation
        #####################
        self.create_workspace(parameters)

        #####################
        # Saving Parameters
        #####################
        self.save_parameters_file(parameters)

        #####################
        # Trajectory Loading
        #####################
        self.timer.start("Trajectory Loading")
        self.trajectoryHandler = TrajectoryHandler(parameters["global"], parameters["matrix"]['parameters'], self.observer)
        self.timer.stop("Trajectory Loading")

        ##############################
        # Distance Matrix Generation
        ##############################
        self.create_matrix(parameters)

        ##############################
        # Actions
        ##############################
        self.perform_actions(parameters)

        self.timer.stop("Global")
        self.notify("Driver Finished", "\n"+str(self.timer))

