"""
Created on 07/09/2012

@author: victor
"""
import random
import numpy
from pyproct.clustering.evaluation.analysis.analysisPopulator import AnalysisPopulator

class BestClusteringSelector(object):

    def __init__(self, parameters):
        """
        Class constructor.

        @param parameters: The global script parameters.
        """
        self.parameters = parameters
        self.criteria = parameters["clustering"]["evaluation"]["evaluation_criteria"]

    def choose_best(self, clustering_info):
        """
        Normalizes the values of the evaluation scores, then calculates the scores for all clusterings and criteria
        and finally chooses the best clustering.

        @param clustering_info: Is the clustering_info structure with clusterings, evaluation info... etc

        @return: The id of the best clustering with the criteria_id with higher score and the score itself.
        """
        if len(clustering_info) == 0:
            print "[WARNING BestClusteringSelector::choose_best] clustering_info is empty."
            return None

        evaluation_types = AnalysisPopulator.get_evaluation_analysis_types(self.parameters)

        # If there were no criteria defined, then the clustering is randomly selected
        if evaluation_types == []:
            return clustering_info[clustering_info.keys()[random.randint(0,len(clustering_info.keys())-1)]]

        for evaluation_type in evaluation_types:
            BestClusteringSelector.normalize_one_evaluation_type(evaluation_type, clustering_info)

        scores = BestClusteringSelector.get_scores_for_all_clusters_and_criterias(self.criteria, clustering_info)


        best_clustering_id, criteria_id, scores = self.get_best_clustering(scores)


        return best_clustering_id, scores

    @classmethod
    def get_best_clustering(cls, scores):
        """
        Selects the clustering with best score.

        @param scores: A scores list as returned by 'get_scores_for_all_clusters_and_criterias'.

        @return: The best clustering Id, the criteria with better results and the score itself.
        """
        best_clustering = (0.0,(None, None))
        for criteria_id in scores:
            for clustering_id in scores[criteria_id]:
                value = scores[criteria_id][clustering_id]
                if value >= best_clustering[0]:
                    best_clustering = (value, (clustering_id,criteria_id))
        best_clustering_id, best_criteria_id = best_clustering[1]
        return best_clustering_id, best_criteria_id, scores

    @classmethod
    def get_scores_for_all_clusters_and_criterias(cls, criteria, clustering_info):
        """
        Calculates all scores for a group of clusterings and criteria.

        @param criteria: Collection of criteria to be applied in the score calculation.

        @param clustering_info: The clustering_info structure with this evaluation_type registered in each of the 'evaluation'
        fields.

        @return: A double dictionary indexed by criteria id and clustering id with all the (clustering id, criteria id) scores.
        """
        scores = {}
        for criteria_id in criteria:
            for clustering_id in clustering_info:
                try:
                    scores[criteria_id][clustering_id] = BestClusteringSelector.get_score_for_criteria(clustering_id,
                                                                                                     clustering_info,
                                                                                                     criteria[criteria_id])
                except KeyError:
                    scores[criteria_id] = {clustering_id : BestClusteringSelector.get_score_for_criteria(clustering_id,
                                                                                                     clustering_info,
                                                                                                     criteria[criteria_id])}
        return scores

    @classmethod
    def get_score_for_criteria(cls, clustering_id, clustering_info, criteria):
        """
        Calculates the score for one clustering and one criteria.

        @param clustering_id: The clustering id of the clustering inside 'clustering_info' we want the score.

        @param clustering_info: The clustering_info structure with this evaluation_type registered in each of the 'evaluation'
        fields.

        @param criteria: Criteria to be applied in the score calculation.

        @return: The score [0..oo).
        """
        evaluation_info = clustering_info[clustering_id]["evaluation"]
        score = 0.0
        accum_weight = 0
        for evaluation_type in criteria:
            value = evaluation_info["Normalized_"+evaluation_type]
            weight = criteria[evaluation_type]["weight"]
            accum_weight += weight
            action = criteria[evaluation_type]["action"]
            if action == ">":
                #Maximize metric
                score += value * weight
            elif action == "<":
                #Minimize metric
                score += (1. - value) * weight
            else:
                print "[ERROR]Criteria action is not valid ( %s )"%action
                exit()
        return score / accum_weight

    @classmethod
    def normalize_one_evaluation_type(cls, evaluation_type, clustering_info):
        """
        Normalizes all the values of one evaluation type in the clustering_info structure in the range [0..1]

        @param evaluation_type: The evaluation type which values we want to normalize.

        @param clustering_info: The clustering_info structure with this evaluation_type registered in each of the 'evaluation'
        fields.

        """
        all_values = []
        for clustering_id in clustering_info:
            all_values.append(clustering_info[clustering_id]["evaluation"][evaluation_type])

        valmax = numpy.max(all_values)
        valmin = numpy.min(all_values)
        for clustering_id in clustering_info:
            value = clustering_info[clustering_id]["evaluation"][evaluation_type]
            if (valmax - valmin) == 0:
                clustering_info[clustering_id]["evaluation"]["Normalized_"+evaluation_type] = 1
            else:
                clustering_info[clustering_id]["evaluation"]["Normalized_"+evaluation_type] =  (value-valmin) / (valmax - valmin)

    @classmethod
    def get_values_for_evaluation_type(cls, evaluation_type, clustering_info):
        """
        Testing helper. Gets all the values of one evaluation type.

        @param evaluation_type: The evaluation type which values we want to recover.

        @param clustering_info: The clustering_info structure with this evaluation_type registered in each of the 'evaluation'
        fields.

        @return: A dictionary indexed by cluster id with the values for the evaluation named 'evaluation_type'
        """
        all_values = {}
        for clustering_id in clustering_info:
            all_values[clustering_id] = clustering_info[clustering_id]["evaluation"][evaluation_type]
        return all_values