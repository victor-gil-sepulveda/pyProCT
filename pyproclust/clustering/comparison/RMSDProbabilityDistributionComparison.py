'''
Created on 20/08/2012

@author: victor
'''
from pyproclust.tools.pdbTools import get_number_of_frames
import numpy
import matplotlib.pyplot as plt
import math

def smoothed(distribution,small_value = 1.0e-8):
    total_number_of_samples = len(distribution)
    samples_in_distrib = numpy.count_nonzero(distribution)
    pc = small_value * (total_number_of_samples - samples_in_distrib) / samples_in_distrib
    smoothed_distrib = numpy.empty(len(distribution))
    for i in range(len(distribution)):
        if distribution[i] == 0:
            smoothed_distrib[i] = small_value
        else:
            smoothed_distrib[i] = distribution[i] - pc
    return smoothed_distrib

def Kullback_Leibler_Divergence(dist1,dist2):
    kl = 0;
    for i in range(len(dist1)):
        try:
            kl += dist1[i] * math.log(dist1[i]/dist2[i],2)
        except Exception:
            print  dist1[i], " ",dist1[i]/dist2[i]
    return kl
    
class RMSDProbabilityDistributionComparison(object):

    def __init__(self, pdb1, pdb2, condensed_distance_matrix, output_file):
        number_of_models_1 = get_number_of_frames(pdb1)
        number_of_models_2 = get_number_of_frames(pdb2)
        max_trajs = max(condensed_distance_matrix.get_data())
        min_trajs = min(condensed_distance_matrix.get_data())
        bin_range = (min_trajs,max_trajs)
        NUM_BINS = 100
        prob_histogram1,bins1 = self.get_probability_histogram(self.get_matrix_data(condensed_distance_matrix,0,number_of_models_1),bin_range,NUM_BINS)
        smoothed_prob_histogram1 = smoothed(prob_histogram1)
        prob_histogram2,bins2 = self.get_probability_histogram(self.get_matrix_data(condensed_distance_matrix,number_of_models_1,number_of_models_2),bin_range,NUM_BINS)
        smoothed_prob_histogram2 = smoothed(prob_histogram2)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pl1 = ax.plot(bins1[:NUM_BINS],smoothed_prob_histogram1, 'b--', linewidth=2)
        pl2 = ax.plot(bins2[:NUM_BINS],smoothed_prob_histogram2, 'r--', linewidth=2)
        ax.grid(True)
        ax.legend([plt.Rectangle((0, 0), 1, 1, fc="b"),plt.Rectangle((0, 0), 1, 1, fc="r")],[pdb1,pdb2])
        plt.savefig(output_file+".png")
        
        self.kl1 = Kullback_Leibler_Divergence(smoothed_prob_histogram1,smoothed_prob_histogram2)
        self.kl2 = Kullback_Leibler_Divergence(smoothed_prob_histogram2,smoothed_prob_histogram1)
        
        open(output_file+".txt","w").write("KL Divergence\nt1 on t2: %f\nt2 on t1: %f"%(self.kl1,self.kl2))
        
    def get_matrix_data(self,matrix,initial_element,number_of_elements):
        traj_data = numpy.empty(number_of_elements*(number_of_elements-1)/2)
        final_element = initial_element+number_of_elements
        m_i = 0
        for i in range(initial_element,final_element-1):
            for j in range(i+1,final_element):
                traj_data[m_i] = matrix[i,j]
                m_i = m_i+1
        return traj_data
    
    def get_probability_histogram(self,data,bin_range,num_bins):
        hist = numpy.histogram(data, num_bins, bin_range)
        float_hist = numpy.asarray(hist[0],dtype=numpy.float32)
        return float_hist/len(data), hist[1]
    
    def getKLnumbers(self):
        return (self.kl1, self.kl2)
    
#    def get_matrix_data(self,matrix):
#        non_repeated_data = [0]*len(matrix.get_data())
#        m_i = 0
#        
#        # Traj1 VS Itself
#        for i in range(self.traj_1_numbr_of_models-1):
#            for j in range(i+1,self.traj_1_numbr_of_models):
#                non_repeated_data[m_i] = matrix[i,j]
#                m_i += 1
#        
#        # Traj2 VS Itself
#        for i in range(self.traj_1_numbr_of_models,self.traj_2_numbr_of_models-1):
#            for j in range(self.traj_1_numbr_of_models+i+1,self.traj_2_numbr_of_models):
#                non_repeated_data[m_i] = matrix[i,j]
#                m_i += 1
#        
#        # Traj1 VS Traj2
#        for i in range(self.traj_1_numbr_of_models):
#            for j in range(self.traj_1_numbr_of_models,self.traj_2_numbr_of_models):
#                non_repeated_data[m_i] = matrix[i,j]
#                m_i += 1
#            
#        return numpy.array(non_repeated_data[:m_i])