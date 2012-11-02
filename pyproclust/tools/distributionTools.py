'''
Created on 16/05/2012

@author: victor
'''
import numpy

def get_distances_for_elems(elems,center,condensed_distance_matrix):
    distances = []
    for e in elems:
        distances.append(condensed_distance_matrix[center,e])
    return distances

def get_distance_std_dev_for_elems(elems,center,condensed_distance_matrix):
    distances = []
    for e in elems:
        distances.append(condensed_distance_matrix[center,e])
    return numpy.std(distances)

def calculate_distribution_area(distribution,bin_size):
    area = 0
    for i in range(len(distribution)):
        area += distribution[i] * bin_size
    return area

def translate_distribution(distribution, bin_distance): 
    new_distr = [0]*len(distribution)
    
    for  i in range(len(distribution)):
        if i-bin_distance>0 and i-bin_distance<len(distribution):
#            print i, i-bin_distance
            new_distr[i] = distribution[i-bin_distance] # to translate we need to sum the opposite
    return new_distr

def calc_distribution_overlap(hist1,hist2,total_bins):
    overlapping = []
    for i in range(total_bins):
        overlapping.append(min(hist1[i],hist2[i]))
    return overlapping

def desp_sign(val1,val2):
    """
    Sign of the step to go from val2 to val1
    """
    if val1 >= val2:
        return 1
    else:
        return -1

def get_nonzero_range(distribution):
    start = 0
    end = len(distribution)-1
    for i in range(len(distribution)):
        if distribution[i]!=0:
            start = i
            break
    
    reverse_range = range(len(distribution))
    reverse_range.reverse()
    for i in reverse_range:
        if distribution[i]!=0:
            end = i
            break
    return start, end



def fit_translation_of_distributions(dist1,dist2,total_bins,bin_size):
    """
    Translates dist2 in order to get the best overlapping (just trying different translations)
    -dist means a displacement to the left
    """
    # Obtain ranges
    start1,end1 = get_nonzero_range(dist1)
#    print start1,end1
    start2,end2 = get_nonzero_range(dist2)
#    print start2,end2
    move_direction   = desp_sign(start1 ,start2)
#    print "move_direction",move_direction
    overlaps = []
    #common range
    search_range = abs(min(start1,start2) - max(end1,end2))
    disp = 0
#    print "s range",search_range
    for i in range(search_range+1):
#        print "disp",disp
        translated = translate_distribution(dist2, disp)
        overlap = calc_distribution_overlap(dist1, translated, total_bins)
        area = calculate_distribution_area(overlap,bin_size)
        overlaps.append((area,disp))
        disp+=move_direction
    del i
#    print overlaps
    return max(overlaps)

