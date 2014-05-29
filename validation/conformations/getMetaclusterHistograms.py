import sys
import json
from pyproct.tools.commonTools import convert_to_utf8
from pyproct.clustering.clustering import Clustering
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt

# def gen_data(cluster, index_to_interpolation_map):
#     data = {}
#     for element in cluster.all_elements:
#         a,b = index_to_interpolation_map[element]
#         if a in data:
#             data[a].append(b)
#         else:
#             data[a] = [b]
#
#     return data

def gen_data(cluster, index_to_interpolation_map,original_frame_number, number_of_interpolated_frames):
    percents = []
    for element in cluster.all_elements:
        try:
            cluster.percents[index_to_interpolation_map[element]] += 1
            a,b = index_to_interpolation_map[element]
            cluster.percents[(b,a)] += 1
        except KeyError:
            cluster.percents[index_to_interpolation_map[element]] = 1
            a,b = index_to_interpolation_map[element]
            cluster.percents[(b,a)] = 1

    for percent in cluster.percents:
        cluster.percents[percent] /= 20.


def get_best_clustering(results_file):
    """
    Loads and returns the best clustering from a results file.
    """
    results = convert_to_utf8(json.loads(open(results_file).read()))
    best_clustering_id =results["best_clustering"]
    best_clustering_dic = results["selected"][best_clustering_id]
    return Clustering.from_dic(best_clustering_dic["clustering"])

def generate_index_interpolation_map(original_frame_number, number_of_interpolated_frames):
    """
    Generates the index-interpolation map i.e. a map that relates dataset element indices with the interpolation
    to which it belongs. It uses the same process that was used to process the frames.
    Ex. index_to_interpolation[0] = (0,3) means that element 0 belongs to the interpolation 0->3

    @param original_frame_number: is the number of frames the original file had (the file from which we calculated the interpolations).
    @param number_of_interpolated_frames: The number of frames we obtain from each interpolation.

    @return: The map.
    """
    acc = 0
    index_to_interpolation = {}
    for i in range(0, original_frame_number-1):
        for j in range(i+1, original_frame_number):
            for k in range(number_of_interpolated_frames):
                index_to_interpolation[acc] = (i,j)
                acc += 1
    print acc
    return index_to_interpolation

if __name__ == '__main__':
    best_clustering = get_best_clustering(sys.argv[1])
    number_of_original_file_frames = int(sys.argv[2])
    number_of_interpolated_frames = int(sys.argv[3])
    index_to_interpolation_map = generate_index_interpolation_map(number_of_original_file_frames, number_of_interpolated_frames)


    for cluster in best_clustering.clusters:
        id = cluster.id
        fig = plt.figure()
        #ax.autoscale(False)
        colors = iter(cm.rainbow(np.linspace(0, 1, number_of_original_file_frames)))
        data = gen_data(cluster, index_to_interpolation_map)

        labels = []
        all_data = []
        for datum in data:
            all_data.append(data[datum])

        labels = [str(d) for d in data.keys()]
        print data
        print all_data
        colors = cm.rainbow(np.linspace(0, 1, number_of_original_file_frames))[0:len(data.keys())]
        print len(data)
        print len(colors)
        counts, bins, patches = plt.hist(all_data, bins = number_of_original_file_frames, range = [0,number_of_original_file_frames],\
                   normed = True, histtype='bar', align='mid', rwidth = 0.50, label = labels)
        #plt.ylim(0, 1. )
        plt.xlabel(id)
        plt.legend()
        plt.xlim(0, number_of_original_file_frames )
        plt.show()
