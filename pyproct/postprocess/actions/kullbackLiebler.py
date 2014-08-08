"""
Created on 16/07/2014

@author: victor
"""
import numpy
import json
import math
import matplotlib.pyplot as plt

class KullbackLieblerPostAction(object):
    KEYWORD = "kullback_liebler"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):
        if KullbackLieblerPostAction.KEYWORD in postprocessing_parameters:
            pass
#             klDiv = KullbackLeiblerDivergence(self.trajectoryHandler.pdbs, self.matrixHandler.distance_matrix)
#             kl_file_path = os.path.join(self.workspaceHandler["matrix"], "kullback_liebler_divergence")
#             klDiv.save(kl_file_path)
#             matrix_image_file_path = os.path.join(self.workspaceHandler["matrix"], parameters["matrix"]["image"]["filename"])
#             self.generatedFiles.append({"description":"Kullback-Leibler divergence",
#                                         "path":matrix_image_file_path,
#                                         "type":"text"})
def smoothed(distribution,small_value = 1.0e-8):
    """
    Applies a smoothing process to the distribution.
    See http://mathoverflow.net/questions/72668/how-to-compute-kl-divergence-when-pmf-contains-0s
    for an explanation about the problem and the solution.
     
    @param distribution: distribution to be smoothed
    @param small_value: value to be set to those bins with 0 probability
     
    @return: The smoothed distribution.
    """
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
 
class KullbackLeiblerDivergence(object):
     
    # Number of bins that the histogram will have when calculating the
    # distribution
    NUM_BINS = 200
     
    def __init__(self, pdb_info, condensedMatrix):
        """
        Class constructor. Does the actual calculation.
         
        @param pdb_info: An structure containing paths and other useful data about the pdbs being used.
        @param condensedMatrix: The actual calculated matrix.
        """
        self.pdb_info = pdb_info
         
        # Getting submatrices.
        submatrices = {}
        next_first_element = 0
        for pdb in pdb_info:
            submatrices[pdb["source"]] = KullbackLeiblerDivergence.get_matrix_data(condensedMatrix, next_first_element, pdb["conformations"])
            next_first_element += pdb["conformations"]
         
        # Getting max and min values
        max_vals = []
        min_vals = []
        for pdb in pdb_info:
            max_vals.append(numpy.max(submatrices[pdb["source"]]))
            min_vals.append(numpy.min(submatrices[pdb["source"]]))
         
        distribution_range = (numpy.min(min_vals), numpy.max(max_vals))
         
        # Generate histograms
        self.histograms = {}
        for pdb in pdb_info:
            prob_histogram, bins =  KullbackLeiblerDivergence.get_probability_histogram(submatrices[pdb["source"]],
                                                                                   distribution_range,
                                                                                   KullbackLeiblerDivergence.NUM_BINS)
            self.histograms[pdb["source"]] =(smoothed(prob_histogram), bins)
         
        # Calculate KL values
        self.KL_matrix = numpy.zeros((len(pdb_info),len(pdb_info)))
        for i in range(len(pdb_info)-1):
            for j in range(i+1,len(pdb_info)):
                self.KL_matrix[i][j] = KullbackLeiblerDivergence.kullback_leibler_divergence_calculation(self.histograms[pdb_info[i]["source"]][0],
                                                                                       self.histograms[pdb_info[j]["source"]][0])
                self.KL_matrix[j][i] = KullbackLeiblerDivergence.kullback_leibler_divergence_calculation(self.histograms[pdb_info[j]["source"]][0],
                                                                                       self.histograms[pdb_info[i]["source"]][0])
         
#         self.pdb1 = pdb_info[0]["source"]
#         self.pdb2 = pdb_info[1]["source"]
#         number_of_models_1 = pdb_info[0]["conformations"]
#         number_of_models_2 = pdb_info[1]["conformations"]
#         
#         first_pdb_submatrix = KullbackLeiblerDivergence.get_matrix_data(condensedMatrix,0,number_of_models_1)
#         second_pdb_submatrix = KullbackLeiblerDivergence.get_matrix_data(condensedMatrix, number_of_models_1, number_of_models_2)
#         
#         max_of_submatrices = max(numpy.max(first_pdb_submatrix),numpy.max(second_pdb_submatrix))
#         min_of_submatrices = min(numpy.min(first_pdb_submatrix),numpy.min(second_pdb_submatrix))
#         distribution_range = (min_of_submatrices, max_of_submatrices)
#         
#         prob_histogram1, self.bins1 = KullbackLeiblerDivergence.get_probability_histogram(first_pdb_submatrix,
#                                                                                           distribution_range,
#                                                                                           KullbackLeiblerDivergence.NUM_BINS)
#         
#         prob_histogram2, self.bins2 = KullbackLeiblerDivergence.get_probability_histogram(second_pdb_submatrix,
#                                                                                           distribution_range,
#                                                                                           KullbackLeiblerDivergence.NUM_BINS)
# 
#         self.smoothed_prob_histogram1 = smoothed(prob_histogram1)
#         self.smoothed_prob_histogram2 = smoothed(prob_histogram2)
#         
#         self.kl1 = KullbackLeiblerDivergence.kullback_leibler_divergence_calculation(self.smoothed_prob_histogram1,self.smoothed_prob_histogram2)
#         self.kl2 = KullbackLeiblerDivergence.kullback_leibler_divergence_calculation(self.smoothed_prob_histogram2,self.smoothed_prob_histogram1)
#     
    def save(self, where):
        """
        Saves a plot of the distributions and the actual values of them.
        @param where: The name of the file without extension (".png" will be appended to the final name). 
        """
        image_path = self.plot_distributions(where)
        return self.to_json(where, image_path)
         
    def plot_distributions(self, where):
        """
        Saves a plot of the distributions.
        @param where: The name of the file without extension (".png" will be appended to the final name). 
        """
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#         ax.plot(self.bins1[:self.NUM_BINS],self.smoothed_prob_histogram1, 'b--', linewidth=2)
#         ax.plot(self.bins2[:self.NUM_BINS],self.smoothed_prob_histogram2, 'r--', linewidth=2)
#         ax.grid(True)
#         ax.legend([plt.Rectangle((0, 0), 1, 1, fc="b"),plt.Rectangle((0, 0), 1, 1, fc="r")],[self.pdb1,self.pdb2])
#         plt.savefig(where+".png")
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for pdb in self.pdb_info:
            smoothed_his, bins  =  self.histograms[pdb["source"]] 
            ax.plot(bins[:self.NUM_BINS],smoothed_his, linewidth = 2, label= os.path.basename(pdb["source"]))
 
        ax.grid(True)
#         ax.legend([plt.Rectangle((0, 0), 1, 1, fc="b"),plt.Rectangle((0, 0), 1, 1, fc="r")],[self.pdb1,self.pdb2])
        plt.legend(prop={'size':6})
        plt.savefig(where+".png")
  
        return where+".png"
     
    def to_json(self, where, image_path):
        """
        Saves the K-L values in a text file containing its json representation.
        @param where: The name of the file without extension (".json" will be appended to the final name).
        @param image_path: Place where we have saved the image with the distribution. 
        """
#         pre_json_dic = {"kl1":self.kl1,"kl2":self.kl2,"image":image_path}
#         for pdb in self.pdb_info:
#         open( where+".json","w").write(json.dumps({"kl1":self.kl1,"kl2":self.kl2,"image":image_path}, indent=4, separators=(',', ': ')))
#         return pre_json_dic
        pre_json = []
        for i in range(len(self.pdb_info)-1):
            for j in range(i+1,len(self.pdb_info)):
                pre_json.append({
                                 "A":self.pdb_info[i],
                                 "B":self.pdb_info[j],
                                 "KL": self.KL_matrix[i][j] 
                                 })
                pre_json.append({
                                 "A":self.pdb_info[j],
                                 "B":self.pdb_info[i],
                                 "KL": self.KL_matrix[j][i] 
                                 })
                 
        open( where+".json","w").write(json.dumps(pre_json, indent=4, separators=(',', ': ')))
     
    @classmethod
    def get_probability_histogram(cls,data,bin_range,num_bins):
        """
        Creates the histogram using numpy.
        @param data: data from which we will create the histogram (in this case the data of one of the submatrices)
        @param bin_range: tuple with the maximum and minimum value the sum of all distributions (not only this one) can have.
        @param num_bins: Number of discrete parts of the distribution.
        """
        hist = numpy.histogram(data, num_bins, bin_range)
        float_hist = numpy.asarray(hist[0],dtype=numpy.float32)
        return float_hist/len(data), hist[1]
     
    @classmethod
    def get_matrix_data(cls,matrix,initial_element,number_of_elements):
        """
        A matrix generated from the concatenation of two pdbs may have 3 submatrices. First is the pairwise matrix of 
        the first pdb, second the pairwise matrix of the second, and third is the distance matrix of pdb1 vs pdb2 
        (in which info is duplicated, as its itself a pairwise distance matrix). This function grabs the data 
        of one of the two first submatrices.
         
        @param matrix: The matrix we are talking about :/
        @param initial_element: is the index of the initial element of the pdb we want to extract the data. For instance
        if we are working with 2 trajectories of 3 and 4 frames, indexes are [tr1:[0,1,2] tr2:[3,4,5,6]], so to extract 
        the first submatrix data, initial_element would be 0 and number_of_elements 3. Extracting the second submatrix will
        need a initial_element value of 3 and a number_of_elements value of 4.
        @param number_of_elements: As explained above, the number of models the pdb we are working with has.
         
        @return: A 1D numpy.array containing the submatrix data. 
        """
        traj_data = numpy.empty(number_of_elements*(number_of_elements-1)/2)
        final_element = initial_element+number_of_elements
        m_i = 0
        for i in range(initial_element, final_element):
            for j in range(i, final_element):
                if(i!=j):
                    traj_data[m_i] = matrix[i,j]
                    m_i = m_i+1
        return traj_data
     
    @classmethod
    def kullback_leibler_divergence_calculation(cls, dist1, dist2):
        """
        Calculates the Kullback - Leibler divergence of two distributions.
        @param dist1: first distribution
        @param dist2: second distribution
         
        @return: The Kullback-Leibler Divergence value
        """
        kl = 0;
        for i in range(len(dist1)):
            try:
                kl += dist1[i] * math.log(dist1[i]/dist2[i],2)
            except ArithmeticError:
                print "dist1[i]", dist1[i],"dist2[i]", dist2[i]
        return kl
     
    def get_calculated_KL_values(self):
        """
        A simple getter...
        """
        return (self.kl1, self.kl2)