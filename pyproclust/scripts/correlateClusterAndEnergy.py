'''
Created on 26/07/2012

@author: victor
'''
import sys
import pickle
import pyproclust.tools.scriptTools as scripts_common

condensed_distance_matrix_path = sys.argv[1]
condensed_distance_matrix = scripts_common.load_matrix(condensed_distance_matrix_path)

energy_file_path = sys.argv[2]
energy_file = open(energy_file_path,'r')

energy_lines = energy_file.readlines()
energy = []
for l in energy_lines:
    energy.append(float(l))
del energy_lines
energy_file.close()

clustering_file_path = sys.argv[3]
clustering_file = open(clustering_file_path,'r')
clustering = pickle.load(clustering_file)
clustering_file.close()

print "Loaded:", clustering.details

structures_to_avoid =  []
i = 0
for e in energy:
    if e > int(sys.argv[4]):
        structures_to_avoid.append(i)
    i = i+1

weighted_mean_position = 0.
mean_energy_difference_to_medoid = 0.
mean_energy_difference_from_medoid = 0.
for cluster in clustering.clusters:
    energy_per_element = []
    center = cluster.calculate_medoid(condensed_distance_matrix)
    for element in cluster.all_elements:
        if not element in structures_to_avoid:
            energy_per_element.append((energy[element],element))
    energy_per_element.sort()
    medoid_index = energy_per_element.index((energy[center],center))
    percent =  medoid_index*100./float(len(cluster.all_elements))
    energy_difference_to_medoid = energy_per_element[medoid_index][0]-energy_per_element[0][0]
    energy_difference_from_medoid = energy_per_element[-1][0] - energy_per_element[medoid_index][0]
    mean_energy_difference_to_medoid += energy_difference_to_medoid
    mean_energy_difference_from_medoid += energy_difference_from_medoid
    weighted_mean_position += len(cluster.all_elements)*percent

mean_energy_difference_to_medoid /=  len(clustering.clusters)   
mean_energy_difference_from_medoid /=  len(clustering.clusters)   
weighted_mean_position /= clustering.total_number_of_elements

print "Mean energy difference to medoid ", mean_energy_difference_to_medoid
print "Mean energy difference from medoid ", mean_energy_difference_from_medoid
print "Weighted position ", weighted_mean_position

## max energy difference
energy.sort()
print "Max energy difference", energy[-100] - energy[0]