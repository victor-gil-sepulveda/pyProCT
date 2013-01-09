'''
Created on 09/01/2013

@author: victor
'''

class CohesionCalculator(object):
    def __init__(self):
        pass
    
    def evaluate(self, cluster, condensed_distance_matrix):
        """
        Returns the cohesion value of a cluster. The condensed matrix given as 
        parameter stores the distances of the elements of the dataset used to extract
        the cluster.
        """
        if cluster.prototype == None:
            return self.__noproto_eval(cluster,condensed_distance_matrix)
        else:
            return self.__proto_eval(cluster,condensed_distance_matrix)
    
    def __noproto_eval(self,cluster,condensed_distance_matrix):
        """
        Evaluation of the cohesion value for a cluster without prototype.
        
        The definition of cohesion would be weight*2*cohesion in the case we follow
        the exact formula in the book []. As we are going to do comparisons, the x2 global
        multiplication doesn't affect.
        
        Cohesion of a cluster of 1 element should be infinite instead of 0...
        """
        size = cluster.get_size()
        weight = 1. / size
        cohesion = 0.
        for i in range(size-1):
            for j in range(i+1,size):
                cohesion = cohesion + condensed_distance_matrix[cluster[i],cluster[j]]
        return weight*cohesion
    
    def __proto_eval(self,cluster,condensed_distance_matrix):
        """
        Evaluation of the cohesion value for a cluster without protorype.
        """
        print "CohesionCalculator,__proto_eval Not implemented"
        exit(-1)
        