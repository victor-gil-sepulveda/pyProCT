"""
Created on 06/06/2012

@author: victor
"""
import numpy
from pyproct.clustering.cluster import get_cluster_sizes
from pyproct.clustering.evaluation.analysis.analysis import Analysis
from pyproct.clustering.evaluation.metrics.graphMetrics import NCut, RatioCut, MinMaxCut
from pyproct.clustering.evaluation.metrics.pcaMetrics import PCAMetric
from pyproct.clustering.evaluation.metrics.cython.normNCut import CythonNCut
from pyproct.clustering.evaluation.metrics.cython.boundedCohesion import CythonMirrorCohesionCalculator
from pyproct.clustering.evaluation.metrics.cython.silhouette import CythonSilhouetteCoefficientCalculator
from pyproct.clustering.evaluation.metrics.cython.meanMinimumDistance import CythonMeanMinimumDistanceCalculator
from pyproct.clustering.evaluation.metrics.boundedCohesion import MirrorCohesionCalculator
from pyproct.clustering.evaluation.metrics.silhouette import SilhouetteCoefficientCalculator
from pyproct.clustering.evaluation.metrics.cohesion import CohesionCalculator
from pyproct.clustering.evaluation.metrics.meanMinimumDistance import MeanMinimumDistanceCalculator
from pyproct.clustering.evaluation.metrics.CalinskiHarabasz import CalinskiHarabaszCalculator
from pyproct.clustering.evaluation.metrics.Dunn import DunnCalculator
from pyproct.clustering.evaluation.metrics.DaviesBouldin import DaviesBouldinCalculator
from pyproct.clustering.evaluation.metrics.gaussianSeparation import GaussianSeparationCalculator
from pyproct.clustering.evaluation.metrics.separation import SeparationCalculator
from pyproct.clustering.evaluation.metrics.compactness import CompactnessCalculator


class AnalysisPopulator(object):

    def __init__(self, matrix_handler, trajectory_handler, parameters):
        """
        Class creator.

        @param matrix_handler: Is the matrix handler object, or any object containing a .distance_matrix value.

        @param trajectory_handler: Is the trajectory handler object, with structural info.

        @param parameters: The script parameters containing a correct 'evaluation' section.
        """
        self.build_all_analysis(matrix_handler, trajectory_handler)
        self.parameters = parameters

    def build_all_analysis(self, matrix_handler, trajectory_handler):
        """
        Create an instance of all possible analysis. If you're not able to create it... you're not able to use it.

        @param matrix_handler: Is the matrix handler object, or any object containing a .distance_matrix value.

        @param trajectory_handler: Is the trajectory handler object, with structural info.
        """
        distance_matrix = matrix_handler.distance_matrix

        self.all_possible_analysis = {}

        # Pure queries
        self.all_possible_analysis["Details"] = Analysis("Details", self.analysis_function_details)
        self.all_possible_analysis["NumClusters"] = Analysis("Number of clusters", self.analysis_function_num_clusters)
        self.all_possible_analysis["NumClusteredElems"] = Analysis("Number of clustered elements", self.analysis_function_total_elements)
        self.all_possible_analysis["MeanClusterSize"] = Analysis("Mean cluster size", self.analysis_function_mean_cluster_size)
        self.all_possible_analysis["PercentInTop4"] = Analysis("Percent in top 4 clusters", self.analysis_function_top_4)
        self.all_possible_analysis["PercentInTop"] = Analysis("Percent in top cluster", self.analysis_function_top_percent)
        self.all_possible_analysis["ClustersTo90"] = Analysis("Clusters to 90", self.analysis_function_num_clusters_to_percent, 90)
        self.all_possible_analysis["NoiseLevel"] = Analysis("Noise level", self.analysis_function_noise_level, distance_matrix.row_length)

        # Evaluators
        self.all_possible_analysis["MirrorCohesion"] = Analysis("MirrorCohesion", self.evaluate_with_calculator,
                                                          {"class":MirrorCohesionCalculator,"matrix":distance_matrix})

        self.all_possible_analysis["Cohesion"] = Analysis("Cohesion", self.evaluate_with_calculator,
                                                          {"class":CohesionCalculator,"matrix":distance_matrix})

        self.all_possible_analysis["Separation"] = Analysis("Separation", self.evaluate_with_calculator,
                                                            {"class":SeparationCalculator,"matrix":distance_matrix})
        self.all_possible_analysis["MinimumMeanSeparation"] = Analysis("MinimumMeanSeparation", self.evaluate_with_calculator,
                                                                       {"class":MeanMinimumDistanceCalculator,"matrix":distance_matrix})
        self.all_possible_analysis["Silhouette"] = Analysis("Silhouette", self.evaluate_with_calculator,
                                                            {"class":SilhouetteCoefficientCalculator,"matrix":distance_matrix})
        self.all_possible_analysis["Calinski-Harabasz"] = Analysis("Calinski-Harabasz", self.evaluate_with_calculator,
                                                                   {"class":CalinskiHarabaszCalculator,"matrix":distance_matrix})
        self.all_possible_analysis["Dunn"] = Analysis("Dunn", self.evaluate_with_calculator,
                                                      {"class":DunnCalculator,"matrix":distance_matrix})
        self.all_possible_analysis["Davies-Bouldin"] = Analysis("Davies-Bouldin", self.evaluate_with_calculator,
                                                                {"class":DaviesBouldinCalculator,"matrix":distance_matrix})
        self.all_possible_analysis["GaussianSeparation"] = Analysis("GaussianSeparation", self.evaluate_with_calculator,
                                                                    {"class":GaussianSeparationCalculator,"matrix":distance_matrix})
        self.all_possible_analysis["Compactness"] = Analysis("Compactness", self.evaluate_with_calculator,
                                                                    {"class":CompactnessCalculator,"matrix":distance_matrix})

        # Cython
        self.all_possible_analysis["CythonMirrorCohesion"] = Analysis("CythonMirrorCohesion", self.evaluate_with_calculator,
                                                                      {"class":CythonMirrorCohesionCalculator,"matrix":distance_matrix})
        self.all_possible_analysis["CythonMinimumMeanSeparation"] = Analysis("CythonMinimumMeanSeparation", self.evaluate_with_calculator,
                                                                             {"class":CythonMeanMinimumDistanceCalculator,"matrix":distance_matrix})
        self.all_possible_analysis["CythonSilhouette"] = Analysis("CythonSilhouette", self.evaluate_with_calculator,
                                                                  {"class":CythonSilhouetteCoefficientCalculator,"matrix":distance_matrix})

        # Graph
        self.all_possible_analysis["RatioCut"] = Analysis("RatioCut", self.evaluate_with_calculator,
                                                          {"class":RatioCut,"matrix":distance_matrix})
        self.all_possible_analysis["NCut"] = Analysis("NCut", self.evaluate_with_calculator,
                                                      {"class":NCut,"matrix":distance_matrix})
        self.all_possible_analysis["NormNCut"] = Analysis("NormNCut", self.analysis_function_norm_n_cut,distance_matrix)
        self.all_possible_analysis["MinMaxCut"] = Analysis("MinMaxCut", self.evaluate_with_calculator,
                                                           {"class":MinMaxCut,"matrix":distance_matrix})

        # Cython & Graph
        self.all_possible_analysis["CythonNormNCut"] = Analysis("CythonNormNCut", self.analysis_function_cython_norm_n_cut,distance_matrix)

        # PCA
        self.all_possible_analysis["PCAanalysis"] = Analysis("PCAanalysis", self.analysis_function_pca, trajectory_handler)

    def get_analysis_list(self):
        """
        Generates the list of required analyzers.
        """
        analysys_list = []

        analysis_types = AnalysisPopulator.get_query_and_evaluation_analysis_types(self.parameters)

        for analysis_type in analysis_types:
            if analysis_type in self.all_possible_analysis:
                analysys_list.append(self.all_possible_analysis[analysis_type])
            else:
                print "[WARNING]", analysis_type, "is not an allowed analysis type"

        return analysys_list

    @classmethod
    def get_query_and_evaluation_analysis_types(self, parameters):
        """
        Returns a list formed by all the analysis we need to calculate, without repetition, from the data in parameters.

        @param parameters: The script parameters.
        """
        queries = parameters["clustering"]["evaluation"]["query_types"]
        queries.extend(AnalysisPopulator.get_evaluation_analysis_types(parameters))
        return list(set(queries))

    @classmethod
    def get_evaluation_analysis_types(self, parameters):
        """
        Returns a list formed by the evaluation types present in criteria.

        @param parameters: The script parameters.
        """
        eval_types =[]
        for evaluation_criteria_id in parameters["clustering"]["evaluation"]["evaluation_criteria"]:
#             for subcriteria in parameters["clustering"]["evaluation"]["evaluation_criteria"][evaluation_criteria_id]:
#                 eval_types.append(subcriteria)
            eval_types.extend(parameters["clustering"]["evaluation"]["evaluation_criteria"][evaluation_criteria_id].keys())
        return list(set(eval_types))

    ################################################################################
    # Next we'll find the analysis functions we have programmed until now.
    # This is the place to make changes if we want to add new analysis functions..
    ################################################################################
    def analysis_function_details(self,clustering):
        """
        Returns the 'details' field of a clustering.
        """
        return clustering.details

    def analysis_function_num_clusters(self,clustering):
        """
        Returns the number of cluster a clustering has.
        """
        return len(clustering.clusters)

    def analysis_function_total_elements(self,clustering):
        """
        Returns the number of elements that are clusterized in this clustering (which
        may not be the total number of elements of the dataset if there were noisy elements)
        """
        return clustering.total_number_of_elements

    def analysis_function_top_4(self,clustering):
        """
        Returns the percentage of elements of the clustering that are in the 4 bigger clusters.
        """
        clustering.sort_clusters_by_size()
        total = 0
        percents = clustering.get_population_percent_of_n_bigger_clusters(4)
        for p in percents:
            total = total+p
        return total

    def analysis_function_num_clusters_to_percent(self,clustering,percent):
        """
        Returns the maximum number of clusters needed to have a percent of the total
        number of clustered elements.
        """
        return clustering.number_of_clusters_to_get_percent(percent)

    def analysis_function_top_percent(self,clustering):
        """
        Returns the percent of elements over the total number of elements of the clustering, that
        have been clustered into the bigger cluster.
        """
        clustering.sort_clusters_by_size()
        return clustering.get_population_percent_of_cluster(0)

    def analysis_function_noise_level(self, clustering, total_elements):
        """
        Returns the percent of noise elements in the dataset.
        """
        return 100.-(clustering.total_number_of_elements/float(total_elements))*100.

    def analysis_function_mean_cluster_size(self,clustering):
        """
        Returns the mean cluster size.
        """
        sizes = get_cluster_sizes(clustering.clusters)[1]
        return numpy.mean(sizes)

    def evaluate_with_calculator(self, clustering, key_args):
        """
        Creates a calculator using a class.
        @param clustering: Clustering to be evaluated.
        @param key_args: A dictionary containing two mandatory keys: 'class' with the class of the calculator to be use and
        'matrix', which contains the distance matrix
        @return: The valua of the evaluation.
        """
        calculator = key_args['class']()
        return calculator.evaluate(clustering, key_args['matrix'])

    def analysis_function_pca(self,clustering, trajectory_handler):
        calculator = PCAMetric(trajectory_handler)
        return calculator.evaluate(clustering)

    def analysis_function_cython_norm_n_cut(self,clustering,condensed_matrix):
        calculator = CythonNCut()
        return calculator.evaluate(clustering, condensed_matrix) / len(clustering.clusters)

    def analysis_function_norm_n_cut(self,clustering,condensed_matrix):
        calculator = NCut()
        return calculator.evaluate(clustering, condensed_matrix) / len(clustering.clusters)
