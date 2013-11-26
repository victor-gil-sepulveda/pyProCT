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
from pyproclust.tools import visualizationTools

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

            if parameters["matrix"]["method"] == "rmsd":
                #Save CA mean squared displacement of best cluster
                #TODO: REFACTORING
                try:
                    displacements_path, CA_mean_square_displacements = visualizationTools.generate_CA_displacements_file(best_clustering,
                                                                                                                         self.trajectoryHandler)

                    self.generatedFiles.append({
                                                "description":"Alpha Carbon mean square displacements",
                                                "path":displacements_path,
                                                "type":"text"
                    })

                    open(displacements_path,"w").write(json.dumps(CA_mean_square_displacements,
                                                          sort_keys=False,
                                                          indent=4,
                                                          separators=(',', ': ')))
                except Exception:
                    print "[ERROR][Driver::postprocess] Impossible to calculate CA displacements file."

            if parameters["matrix"]["method"] == "distance":
                try:
                    centers_path, centers_contents = visualizationTools.generate_selection_centers_file(parameters,
                                                                                                        best_clustering,
                                                                                                        self.workspaceHandler,
                                                                                                        self.trajectoryHandler)

                    self.generatedFiles.append({
                                                "description":"Centers of the selection used to calculate distances",
                                                "path":centers_path,
                                                "type":"text"
                    })

                    open(centers_path,"w").write(json.dumps(centers_contents,
                                              sort_keys=False,
                                              indent=4,
                                              separators=(',', ': ')))
                except Exception:
                    print "[ERROR][Driver::postprocess] Impossible to calculate selection centers file."


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

