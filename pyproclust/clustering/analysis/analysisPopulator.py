'''
Created on 06/06/2012

@author: victor
'''
import numpy
from pyproclust.clustering.cluster import get_cluster_sizes
from pyproclust.clustering.analysis.analysis import Analysis
from pyproclust.clustering.metrics.cython.CythonMetrics import CythonBoundedCohesionCalculator,\
    CythonSilhouetteCoefficientCalculator,CythonMeanMinimumDistanceCalculator
from pyproclust.clustering.metrics.graphMetrics import NCut, RatioCut, MinMaxCut
from pyproclust.clustering.metrics.pcaMetrics import PCAMetric


class AnalysisPopulator(object):
    
    def __init__(self,analyzer,distance_matrix,pdb_struct,analysis_type_list):
        ##############################################
        # Create an instance of all possible analysis
        # If you're not able to create it... you're 
        # not able to use it. 
        ##############################################
        self.all_possible_analysis = {}
        
        self.all_possible_analysis["Details"]=Analysis("Details",self.analysis_function_details)
        self.all_possible_analysis["NumClusters"]=Analysis("Number of clusters",self.analysis_function_num_clusters)
        self.all_possible_analysis["NumClusteredElems"]=Analysis("Number of clustered elements",self.analysis_function_total_elements)
        self.all_possible_analysis["MeanClusterSize"]=Analysis("Mean cluster size",self.analysis_function_mean_cluster_size)
        self.all_possible_analysis["PercentInTop4"]=Analysis("Percent in top 4 clusters",self.analysis_function_top_4)
        self.all_possible_analysis["PercentInTop"]=Analysis("Percent in top cluster",self.analysis_function_top_percent)
        self.all_possible_analysis["ClustersTo90"]=Analysis("Clusters to 90",self.analysis_function_num_clusters_to_percent,90)
        if distance_matrix:
            self.all_possible_analysis["Cohesion-Separation"]=Analysis("Cohesion and Separation",self.analysis_function_get_cohesion_separation,[distance_matrix,True])
            self.all_possible_analysis["NoiseLevel"]=Analysis("Noise level",self.analysis_function_noise_level,distance_matrix.row_length)
            self.all_possible_analysis["MirrorCohesion"]=Analysis("MirrorCohesion",self.analysis_function_mirror_bounded_cohesion,distance_matrix)
            self.all_possible_analysis["MinimumMeanSeparation"]=Analysis("MinimumMeanSeparation",self.analysis_function_minimum_mean_distance,distance_matrix)
            self.all_possible_analysis["Silhouette"]=Analysis("Silhouette",self.analysis_function_get_silhouette,distance_matrix)
            # Cython
            self.all_possible_analysis["CythonMirrorCohesion"]=Analysis("CythonMirrorCohesion",self.analysis_function_cython_cohesion,distance_matrix)
            self.all_possible_analysis["CythonMinimumMeanSeparation"]=Analysis("CythonMinimumMeanSeparation",self.analysis_function_cython_mean_distance,distance_matrix)
            self.all_possible_analysis["CythonSilhouette"]=Analysis("CythonSilhouette",self.analysis_function_cython_silhouette,distance_matrix)
            # Graph
            self.all_possible_analysis["RatioCut"]=Analysis("RatioCut",self.analysis_function_ratio_cut,distance_matrix)
            self.all_possible_analysis["NCut"]=Analysis("NCut",self.analysis_function_n_cut,distance_matrix)
            self.all_possible_analysis["NormNCut"]=Analysis("NormNCut",self.analysis_function_norm_n_cut,distance_matrix)
            self.all_possible_analysis["MixMaxCut"]=Analysis("MixMaxCut",self.analysis_function_minmax_cut,distance_matrix)
            # PCA
            self.all_possible_analysis["PCAanalysis"]=Analysis("PCAanalysis",self.analysis_function_pca,pdb_struct)
            
        ##############################################
        # We will add only those required by the parameters
        ##############################################
        analyzer.add_analysis(self.all_possible_analysis["Details"])
        
        for analysis_type in analysis_type_list:
            if analysis_type in self.all_possible_analysis:
                analyzer.add_analysis(self.all_possible_analysis[analysis_type])
            else:
                print "[WARNING]", analysis_type, "is not an allowed analysis type"
   
    ########################################################
    # Next we'll find the analysis functions we have programmed until now
    # This is the place to make changes if we want to add new analysis functions..
    ########################################################
    
    def analysis_function_ratio_cut(self,clustering,condensed_matrix):
        calculator = RatioCut()
        return calculator.evaluate(clustering, condensed_matrix)
    
    def analysis_function_n_cut(self,clustering,condensed_matrix):
        calculator = NCut()
        return calculator.evaluate(clustering, condensed_matrix)
    
    def analysis_function_norm_n_cut(self,clustering,condensed_matrix):
        calculator = NCut()
        return calculator.evaluate(clustering, condensed_matrix) / len(clustering.clusters)
    
    def analysis_function_minmax_cut(self,clustering,condensed_matrix):
        calculator = MinMaxCut()
        return calculator.evaluate(clustering, condensed_matrix)
        
    def analysis_function_pca(self,clustering,pdb_struct):
        calculator = PCAMetric()
        return calculator.evaluate(clustering, pdb_struct)
        
    def analysis_function_cython_cohesion(self,clustering,distance_matrix):
        """
        Analysis using a cython faster version.
        """
        calculator = CythonBoundedCohesionCalculator()
        return calculator.evaluate(clustering,distance_matrix)
    
    def analysis_function_cython_mean_distance(self,clustering,distance_matrix):
        """
        Analysis using a cython faster version.
        """
        calculator = CythonMeanMinimumDistanceCalculator()
        return calculator.evaluate(clustering,30,distance_matrix)
    
    def analysis_function_cython_silhouette(self,clustering,distance_matrix):
        """
        Analysis using a cython faster version.
        """
        calculator = CythonSilhouetteCoefficientCalculator()
        return calculator.evaluate(clustering,distance_matrix)
    
    def analysis_function_num_clusters(self,clusterization):
        """
        Returns the number of cluster a clustering has.
        """
        return len(clusterization.clusters)
    
    def analysis_function_total_elements(self,clusterization):
        """
        Returns the number of elements that are clusterized in this clustering (which
        may not be the total number of elements of the dataset if there were noisy elements)
        """
        return clusterization.total_number_of_elements
    
    def analysis_function_top_4(self,clusterization):
        """
        Returns the percentage of elements of the clustering that are in the 4 bigger clusters.
        """
        clusterization.sort_clusters_by_size()
        total = 0
        percents = clusterization.get_population_percent_of_n_bigger_clusters(4)
        for p in percents:
            total = total+p
        return total
    
    def analysis_function_top_percent(self,clusterization):
        """
        Returns the percent of elements over the total number of elements of the clustering, that 
        have been clustered into the bigger cluster.
        """
        clusterization.sort_clusters_by_size()
        return clusterization.get_population_percent_of_cluster(0)
    
    def analysis_function_num_clusters_to_percent(self,clusterization,percent):
        """
        Returns the maximum number of clusters needed to have a percent of the total 
        number of clustered elements.
        """
        return clusterization.number_of_clusters_to_get_percent(percent)
    
    def analysis_function_get_cohesion_separation(self,clusterization,distance_matrix_unprototype):
        """
        Calculates cohesion and separation of a clustering.
        """
        cohe,sep = clusterization.cohesion_and_separation_factor(distance_matrix_unprototype[0],distance_matrix_unprototype[1])
        return (cohe,sep)
    
    def analysis_function_get_silhouette(self,clusterization,distance_matrix):
        """
        Calculates silhouette factor of a clustering.
        """
        return clusterization.silhouette_factor(distance_matrix)
    
    def analysis_function_details(self,clusterization):
        """
        Returns the 'details' field of a clustering.
        """
        return clusterization.details
    
    def analysis_function_mirror_bounded_cohesion(self,clusterization,distance_matrix):
        """
        Calculates the mirror bounded cohesion of a clustering.
        """
        return 1-clusterization.bounded_cohesion(distance_matrix)
    
    def analysis_function_minimum_mean_distance(self,clusterization,distance_matrix):
        """
        Calculates the minimum mena distance...
        """
        return clusterization.minimum_mean_distance(0.1,distance_matrix) # 10% subsampling
    
    def analysis_function_noise_level(self,clusterization,total_elements):
        """
        Returns the percent of noise elements in the dataset.
        """
        return 100-(clusterization.total_number_of_elements/float(total_elements))*100
    
    def analysis_function_mean_cluster_size(self,clusterization):
        """
        Returns the mean cluster size.
        """
        sizes = get_cluster_sizes(clusterization.clusters)[1]
        return numpy.mean(sizes)