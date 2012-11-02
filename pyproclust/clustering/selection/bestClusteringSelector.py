'''
Created on 07/09/2012

@author: victor
'''

class BestClusteringSelector(object):
    
    def __init__(self, value_maps):
        self.cluster_score_value_maps = value_maps
    
    def calcMaxScoring(self,value_map):    
        """
        Calculates the maximum possible scoring to be used in scoring normalization
        """
        max_scoring = 0.
        for key in value_map:
            max_scoring += value_map[key][0]
        return max_scoring
        
    def chooseBestClustering(self,results_pack):
        """
        Calculates the best scoring (and cluster) given a list off value_maps.
        """
        if len(self.cluster_score_value_maps) == 0:
            print "Error choosing the best clustering: there were no params for choosing. Exiting..."
            exit()
        best_scores = []
        all_scores = []
        for value_map in self.cluster_score_value_maps:
            scores = self._choose_best_clustering(results_pack,value_map)
            all_scores.append((value_map,scores))
            best_scores.append(scores[0])
            print "Best score for ", value_map, " is ", best_scores[-1][0]
        best_scores.sort(reverse=True)
        return best_scores[0], all_scores
        
    def _choose_best_clustering(self,results_pack,value_map):
        """
        Returns the best clustering using the data in the parameters and the results_pack.
        'results_pack' is a list of tuples where each tuple first element is a clustering,
        and the second element is a dictionary indexed by analysis name with the numerical 
        value of this analysis for the given clustering. 
        A score value of 1 says that the best cluster had also the best partial scoring values.
        """
        clustering_scores = []
        for pack in results_pack:
            clustering = pack[0]
            eval_list = pack[1]
            norm_score = self.scoreClustering(eval_list,value_map)
            clustering_scores.append((norm_score,clustering))
        
        # Return the higher one
        clustering_scores.sort(reverse=True)
        return clustering_scores
    
    def scoreClustering(self,eval_list,value_map):
        """
        Calculates the score for a given evaluation list. An evaluation list is a dictionary that 
        represents the score of each of the possible metrics for a clustering (and is a value in [0,1]).
        This functions calculates a [0,1] normalized scoring using this eval list.
        """
        score = 0.
        for value_key in value_map:
            # It the keys are not defined the program has to crash
            pair = value_map[value_key]
            # Example: 
            # { "NCut":(1.0,">") }
            multiplier = pair[0]
            min_or_max = pair[1]
            val = eval_list[value_key]
            #PRECONDITION: Values are normalized to 1
            if min_or_max == ">": # metric to Maximize
                score += val * multiplier
            else: # metric to Minimize
                score += (1-val) * multiplier
        final_score = score/self.calcMaxScoring(value_map)
        return final_score