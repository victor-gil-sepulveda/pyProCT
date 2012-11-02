'''
Created on 12/06/2012

@author: victor
'''
import numpy
import pyproclust.tools.commonTools as common
from pyproclust.tools import distributionTools
from pyproclust.clustering.cluster import Cluster
from pyproclust.tools.pdbTools import get_number_of_frames
from pyproclust.tools.plotTools import ballGraph, plotDataCards, pieChartCreation,\
    plotSummaryTable, barGraphCreation
import Image

def calculate_centers_difference(cluster_type, cluster, condensed_distance_matrix, A_elems = [], B_elems =[] ): #type = 'M'(Mixed) or 'P' (Pure)
    if cluster_type == 'A' or cluster_type == 'B' or cluster_type == 'P':
        return 0;
    else:
        if cluster_type == 'M':
            A_medoid = cluster.calculate_biased_medoid(condensed_distance_matrix, A_elems)
            B_medoid = cluster.calculate_biased_medoid(condensed_distance_matrix, B_elems)
            centers_difference = condensed_distance_matrix[A_medoid,B_medoid] 
            return centers_difference;
        else:
            raise TypeError("Type has to be 'M' or 'P'")

def calculate_cluster_mean_distance_and_dispersion(cluster, selected_elems, condensed_distance_matrix):
    medoid = cluster.calculate_biased_medoid(condensed_distance_matrix, selected_elems)
    # Remove the medoid from the selected_elements
    s_copy = list(selected_elems)
    s_copy.remove(medoid)
    
    if s_copy == []:
        mean_distance = 0
        std_distance = 0
    else:
        distances = distributionTools.get_distances_for_elems(s_copy, medoid, condensed_distance_matrix) #@UnusedVariable
        mean_distance = numpy.mean(distances)
        std_distance = numpy.std(distances)
    return mean_distance, std_distance

def calculate_max_distance(cluster,condensed_distance_matrix):
    medoid = cluster.calculate_biased_medoid(condensed_distance_matrix, cluster.all_elements)
    distances = distributionTools.get_distances_for_elems(cluster.all_elements, medoid, condensed_distance_matrix)
    return numpy.max(distances)

def dic2ListWoZeros(mydic):
    mylist = []
    for key in mydic:
        elem = mydic[key]
        if elem!=0:
            mylist.append(elem)
    return mylist

def getTotalNumberOfElements(list_of_cluster_lists):
    total = 0
    for cluster_list in list_of_cluster_lists:
        for c in cluster_list:      
            total += c.get_size()
    return total

class Separator(object):
    
    def __init__(self):
        pass
    
    def c_separate(self, clusters, traj1_pdb, traj2_pdb):
        # Retrieve number of models. traj1 will be the first trajectory used to get the clustering, and traj2 the second
        # one.
        self.traj_1_numbr_of_models = get_number_of_frames(traj1_pdb)
        self.traj_2_numbr_of_models = get_number_of_frames(traj2_pdb)
        
        # Calculate data ranges inside the clustering for each trajectory
        traj_1_elements,traj_2_elements = common.gen_consecutive_ranges(self.traj_1_numbr_of_models, self.traj_2_numbr_of_models)
        
        # Now separate clusters
        return self.separate_pure_from_mixed_clusters(clusters,traj_1_elements,traj_2_elements) 
    
    def separate(self, clustering, traj1_pdb, traj2_pdb):
        return self.c_separate(clustering.clusters, traj1_pdb, traj2_pdb) 
    
    def separate_pure_from_mixed_clusters(self, clusters, range_1, range_2):
        """
        Separates all the clusters in clutering in 3 categories. The first one would be the 
        clusters which have only elements in range1, the second would be the clusters with only
        elements in range2, and then the third and last will be the formed by the cluster which
        have elements in both ranges at the same time. The mixed clusters list will be a 
        list containing tuples with (cluster_instance, elements of the cluster in range 1, elements in range 2)
        """
        pure_clusters_traj1 = []
        pure_clusters_traj2 = []
        mixed_clusters = []
        traj_1_elements_set = set(range_1)
        traj_2_elements_set = set(range_2)
        for c in clusters:
            elements_set = set(c.all_elements)
            elems_of_c_in_traj1 = elements_set.intersection(traj_1_elements_set)
            elems_of_c_in_traj2 = elements_set.intersection(traj_2_elements_set)
            if len(elems_of_c_in_traj1) != 0 and len(elems_of_c_in_traj2) != 0:
                mixed_clusters.append((c,list(elems_of_c_in_traj1),list(elems_of_c_in_traj2)))
            elif len(elems_of_c_in_traj1) != 0:
                pure_clusters_traj1.append(c)
            else:
                pure_clusters_traj2.append(c)
        return pure_clusters_traj1 , pure_clusters_traj2, mixed_clusters    

class ClusterStatisticalData(object):
    def __init__(self, cluster, cluster_type, condensed_distance_matrix,A_elements = [], B_elements = []):
        self.cluster_type = cluster_type
        self.number_of_elements = cluster.get_size()
        self.center_difference = 0
        self.mean_dist_to_center, self.dispersion = calculate_cluster_mean_distance_and_dispersion(cluster,cluster.all_elements,condensed_distance_matrix)
        self.max_distance = calculate_max_distance(cluster,condensed_distance_matrix)
        
        # If mixed ...
        self.A_part_cluster_data = None
        self.B_part_cluster_data = None
        if cluster_type == 'M':
            self.A_part_cluster = Cluster(None,A_elements)
            self.B_part_cluster = Cluster(None,B_elements)
            self.A_part_cluster_data = ClusterStatisticalData(self.A_part_cluster,'P',condensed_distance_matrix)
            self.B_part_cluster_data = ClusterStatisticalData(self.B_part_cluster,'P',condensed_distance_matrix)
            self.center_difference = calculate_centers_difference(cluster_type, cluster, condensed_distance_matrix,A_elements,B_elements)
    
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        res = ""
        if self.cluster_type == 'A' or self.cluster_type == 'B' or self.cluster_type == 'P':
            res = "[Cluster Stats. Data -\n\tType: 'Pure'\n\tSize: %d\n\tMean Dist. to centre: %.3f\n\tDispersion: %.3f]\n"%(self.number_of_elements,\
                                                                                                                     self.mean_dist_to_center,\
                                                                                                                     self.dispersion)
        else:
            res = "[Cluster Stats. Data -\n\ttype: 'Mixed'\n\tSize: %d\n\tMean Dist. to centre: %.3f\n\tDispersion: %.3f\n\tDist. between centers: %.3f\n\n\tA elements inside mixed:"%(self.number_of_elements,\
                                                                                                                     self.mean_dist_to_center,\
                                                                                                                     self.dispersion,self.center_difference)\
            +str(self.A_part_cluster_data )+"\n\tB elements inside mixed:"+str(self.B_part_cluster_data )+"]\n"
        return res
            
class ClusteringStatisticalAnalyzer(object):
    def __init__(self, clustering, traj_A_pdb, traj_B_pdb, condensed_matrix, trajectory_comparison = True):
        self.traj_A_pdb = traj_A_pdb
        self.traj_B_pdb = traj_B_pdb
        self.matrix = condensed_matrix
        self.clustering = clustering
        self.total_number_of_elements = clustering.total_number_of_elements
        self.trajectory_comparison = trajectory_comparison
        
        if not trajectory_comparison:
            self.pure_A = clustering.clusters
            self.pure_B = []
            self.mixed_clusters_with_elements = []
        else:
            separator = Separator()
            self.pure_A , self.pure_B, self.mixed_clusters_with_elements = separator.separate(clustering, traj_A_pdb, traj_B_pdb)
            self.traj_A_number_of_models = separator.traj_1_numbr_of_models
            self.traj_B_number_of_models = separator.traj_2_numbr_of_models
        
        self.trajectory_comparison = trajectory_comparison
        
        self.cluster_statistical_data = []
        self.statistics_dic = {}
        
        self.mixed_clusters_without_elements = []
        for m in self.mixed_clusters_with_elements:
            self.mixed_clusters_without_elements.append(m[0])
    
    def per_cluster_analytics(self):
        for cluster in self.pure_A:
            self.cluster_statistical_data.append(ClusterStatisticalData(cluster,'A',self.matrix))
        
        for cluster in self.pure_B:
            self.cluster_statistical_data.append(ClusterStatisticalData(cluster,'B',self.matrix))
            
        for cluster, a_elems, b_elems in self.mixed_clusters_with_elements:
            self.cluster_statistical_data.append(ClusterStatisticalData(cluster,'M',self.matrix,a_elems,b_elems))
    
    def per_clustering_analytics(self):
        if self.trajectory_comparison:
            self.statistics_dic["number_elements_pure_A"] = getTotalNumberOfElements([self.pure_A])
            self.statistics_dic["elems_percent_pure_A"]= self.statistics_dic["number_elements_pure_A"]*100./self.total_number_of_elements
            self.statistics_dic["elems_self_percent_pure_A"]= self.statistics_dic["number_elements_pure_A"]*100./self.traj_A_number_of_models
            self.statistics_dic["number_elements_pure_B"] = getTotalNumberOfElements([self.pure_B])
            self.statistics_dic["elems_percent_pure_B"]= self.statistics_dic["number_elements_pure_B"]*100./self.total_number_of_elements
            self.statistics_dic["elems_self_percent_pure_B"]= self.statistics_dic["number_elements_pure_B"]*100./self.traj_B_number_of_models
            self.statistics_dic["number_elements_mixed"] = getTotalNumberOfElements([self.mixed_clusters_without_elements])
            self.statistics_dic["elems_percent_mixed"]= self.statistics_dic["number_elements_mixed"]*100./self.total_number_of_elements
    
class ClusteringFinalReportGenerator(object):
    def __init__(self):
        pass
           
def getAllOfAnAttribute(mylist, attr_string):
    gathering_list = []
    item = None
    for this_object in mylist: #@UnusedVariable
        try:
            item = eval("this_object."+attr_string)
        except:
            item = 0
        gathering_list.append(item)
    return gathering_list

class ClusteringPlotsGenerator(object):
    def __init__(self, clusteringStatisticalAnalyzer):
        self.stats_analyzer = clusteringStatisticalAnalyzer
        self.statistical_data= self.stats_analyzer.cluster_statistical_data
        self.stats_analyzer.cluster_statistical_data.sort(reverse=True, key=lambda csd: csd.number_of_elements)
        self.data_lists = {}
        self.data_lists["cluster_mean_distances"] = ("medi:",getAllOfAnAttribute(self.statistical_data,'mean_dist_to_center'))
        self.data_lists["cluster_max_distances"] =  ("mxdi:",getAllOfAnAttribute(self.statistical_data,'max_distance'))
        self.data_lists["cluster_dispersions"] =    ("disp:",getAllOfAnAttribute(self.statistical_data,'dispersion'))
        self.data_lists["center_differences"]  =    ("cdif:",getAllOfAnAttribute(self.statistical_data,'center_difference'))
        self.data_lists["cluster_sizes"] =          ("size:",getAllOfAnAttribute(self.statistical_data,'number_of_elements'))
        self.cluster_types = getAllOfAnAttribute(self.statistical_data,'cluster_type')
        
        A_data = getAllOfAnAttribute(self.statistical_data,'A_part_cluster_data')
        self.data_lists["A_sizes"] = ("Aele:",getAllOfAnAttribute(A_data,'number_of_elements'))
        B_data = getAllOfAnAttribute(self.statistical_data,'B_part_cluster_data')
        self.data_lists["B_sizes"] = ("Bele:",getAllOfAnAttribute(B_data,'number_of_elements'))
        
        self.str_colors = {'A':'#468C54','B':'#4786A1','M':'#949E5D'}
        self.colors = {'A':(70, 140, 84),'B':(71, 134, 161),'M':(148, 158, 93)}
        
    def generate_and_compose_big_plot(self,composition_size, max_radius, ball_horizontal_separation, ball_vertical_separation):
        pie_chart_percents =     (0.4,0.6)
        ball_graph_percents =    (0.65,1)
        summary_table_percents = (0.6,0.4)
        
        canvas = Image.new("RGBA", composition_size, color=(255,)*4)
        
        ball_plot_size = (int(composition_size[0]*ball_graph_percents[0]),int(composition_size[1]))
        ball_plot =  self.generate_ball_plot(ball_plot_size, max_radius, ball_horizontal_separation, ball_vertical_separation,self.colors)
        graph_plot = None
        summary_table_size = None
        
        if self.stats_analyzer.trajectory_comparison == True:
            fractions = (self.stats_analyzer.statistics_dic["number_elements_pure_A"],\
                         self.stats_analyzer.statistics_dic["number_elements_pure_B"],\
                         self.stats_analyzer.statistics_dic["number_elements_mixed"])
            graph_size = (int(composition_size[0]*pie_chart_percents[0]),int(composition_size[1]*pie_chart_percents[1]))
            graph_plot = pieChartCreation(graph_size, fractions, self.stats_analyzer.traj_A_pdb, self.stats_analyzer.traj_B_pdb,self.str_colors)
            summary_table_size = (int(composition_size[0]*summary_table_percents[0]),int(composition_size[1]*summary_table_percents[1]))
        else:
            summary_table_size = (int(composition_size[0]*summary_table_percents[0]),int(composition_size[1]))
        
        summary_table = plotSummaryTable(summary_table_size, self.stats_analyzer.statistics_dic,\
                                        self.data_lists,\
                                        self.stats_analyzer.total_number_of_elements,\
                                        self.stats_analyzer.trajectory_comparison)
        
        # Compose images
        if self.stats_analyzer.trajectory_comparison == True:
            canvas.paste(graph_plot, (0, 0 )) 
            canvas.paste(summary_table, (0,int(composition_size[1]*pie_chart_percents[1]))) 
            canvas.paste(ball_plot, (int(composition_size[0]*pie_chart_percents[0]),0)) 
        else:
            canvas.paste(summary_table, (0,0)) 
            canvas.paste(ball_plot, (int(composition_size[0]*0.4),0)) 
        
        return canvas 
    
    def generate_and_compose_small_plot(self, composition_size):
        canvas = Image.new("RGBA", composition_size, color=(255,)*4)
        
        if self.stats_analyzer.trajectory_comparison == True:
            pie_chart_percents =      (1,0.7) # Overlapping with the next...
            summary_table_percents=   (1,0.28)
            graph_bar_percents =      (1,0.22)
        else:
            summary_table_percents  = (1,0.5)
            graph_bar_percents =      (1,0.5)
        
        if self.stats_analyzer.trajectory_comparison == True:
            fractions = (self.stats_analyzer.statistics_dic["number_elements_pure_A"],\
                         self.stats_analyzer.statistics_dic["number_elements_pure_B"],\
                         self.stats_analyzer.statistics_dic["number_elements_mixed"])
            graph_size = (int(composition_size[0]*pie_chart_percents[0]),int(composition_size[1]*pie_chart_percents[1]))
            graph_plot = pieChartCreation(graph_size, fractions, self.stats_analyzer.traj_A_pdb, self.stats_analyzer.traj_B_pdb,self.str_colors)
        
        summary_table_size = (int(composition_size[0]*summary_table_percents[0]),\
                                  int(composition_size[1]*(summary_table_percents[1])))
        
        summary_table = plotSummaryTable(summary_table_size, self.stats_analyzer.statistics_dic,\
                                        self.data_lists,\
                                        self.stats_analyzer.total_number_of_elements,\
                                        self.stats_analyzer.trajectory_comparison)
        if self.stats_analyzer.trajectory_comparison == True:
            bar_plot = barGraphCreation(self.data_lists["A_sizes"][1],\
                                        self.data_lists["B_sizes"][1],\
                                        self.data_lists["cluster_sizes"][1],\
                                        self.cluster_types,\
                                        numpy.sum(self.data_lists["cluster_sizes"][1]),self.str_colors,\
                                        (int(composition_size[0]*graph_bar_percents[0]),\
                                         int(composition_size[1]*graph_bar_percents[1])))
        else:
            bar_plot = barGraphCreation(self.data_lists["cluster_sizes"][1],\
                                        self.data_lists["B_sizes"][1],\
                                        self.data_lists["cluster_sizes"][1],\
                                        self.cluster_types,\
                                        numpy.sum(self.data_lists["cluster_sizes"][1]),self.str_colors,\
                                        (int(composition_size[0]*graph_bar_percents[0]),\
                                         int(composition_size[1]*graph_bar_percents[1])))
            
        # Compose images
        if self.stats_analyzer.trajectory_comparison == True:
            canvas.paste(graph_plot, (0, 0 )) 
            canvas.paste(summary_table, (0,int(composition_size[1]*(1-(summary_table_percents[1]+graph_bar_percents[1]))))) 
            canvas.paste(bar_plot, (0,int(composition_size[1]*(1-graph_bar_percents[1])))) 
        else:
            canvas.paste(summary_table, (0,0)) 
            canvas.paste(bar_plot, (0,int(composition_size[1]*summary_table_percents[1]))) 
#        
        return canvas 
    
    def generate_ball_plot(self, ball_plot_size, max_radius, ball_horizontal_separation, ball_vertical_separation,colors):
        # Sort by size
        image_with_balls, real_max_radius = ballGraph(len(self.stats_analyzer.cluster_statistical_data), ball_plot_size,\
                  max_radius, ball_horizontal_separation,\
                  ball_vertical_separation, self.data_lists["cluster_sizes"][1], self.data_lists["cluster_dispersions"][1],\
                  colors,self.cluster_types, self.data_lists["A_sizes"][1],self.data_lists["B_sizes"][1])
        return plotDataCards(image_with_balls, ball_plot_size, len(self.stats_analyzer.cluster_statistical_data),\
                             real_max_radius, ball_horizontal_separation, ball_vertical_separation, self.data_lists,\
                             key_exceptions = ["A_sizes","B_sizes"])
    
    def generate_transition_graph(self,is_a_trajectory_comparison):
        pass
#            ATENCION!!!! ESTO REQUIERE LAS STDDEV DE A CON RESPECTO AL CENTRO DEL CLUSTER!!
#            purge_mixed_clusters_and_do_graph(mixed_clusters_with_elements, pure_A,condensed_distance_matrix, scalc.std_devs_from_A, state_graph_path)
#  
    
        