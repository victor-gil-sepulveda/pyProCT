'''
Created on 18/04/2012

@author: victor
'''
from pyproclust.tools.scriptTools import load_matrix
import sys
import numpy
import matplotlib.pyplot as plt
from pyproclust.matrix.condensedMatrix import CondensedDistanceMatrix

if __name__ == '__main__':
    matrix_file = sys.argv[1]
    matrix_image_file = sys.argv[2]
    condensed_distance_matrix = load_matrix(matrix_file)
    # Normalize
    contents = condensed_distance_matrix.get_data()
    _max = numpy.max(contents)
    _min = numpy.min(contents)
    
    norm_contents = (contents - _min) / (_max - _min)
    del condensed_distance_matrix
    
    norm_condensed = CondensedDistanceMatrix(norm_contents)
    _max = numpy.max(norm_contents)
    _min = numpy.min(norm_contents)
    print _max," ",_min
    
    complete = numpy.zeros([norm_condensed.row_length,norm_condensed.row_length] ,dtype=numpy.float)
    
    for i in range(norm_condensed.row_length-1):
        for j in range(i+1,norm_condensed.row_length):
            complete[i][j] = norm_condensed[i,j]
            complete[j][i] = norm_condensed[i,j]
    

    fig = plt.figure()
    plt.gray()
    ax = fig.add_subplot(111)
    cax = ax.imshow(complete, interpolation='nearest')
#    plt.show()
    plt.savefig(matrix_image_file)