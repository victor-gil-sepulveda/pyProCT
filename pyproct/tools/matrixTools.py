"""
Created on 13/11/2013

@author: victor
"""
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix

def get_submatrix( old_matrix, elements):
    """
    Gets the distances of the elements in 'elements' to build a new distance matrix.
    
    @param old_matrix: The initial matrix.
    @param elements: An array of elements of old_matrix (old_matrix.row_length >= len(elements), 
    max(elements)>=old_matrix.row_length)
    
    @return: The new (probably smaller) matrix.
    """
#     This is not working. TODO: Review getting data from matrix once is changed
#     N = len(elements)
#     new_matrix = CondensedMatrix([1.]*((N*(N-1)/2)))
#     for i in range(len(elements)-1):
#         e_i = elements[i]
#         for j in range(i+1,len(elements)):
#             e_j = elements[j]
#             print i, j
#             new_matrix[i, j] = old_matrix[e_i, e_j]
#     return new_matrix
    N = len(elements)
    inner_data = numpy.zeros((N*(N-1))/2)
    k = 0
    for i in range(len(elements)-1):
        e_i = elements[i]
        for j in range(i+1,len(elements)):
            e_j = elements[j]
            inner_data[k] = old_matrix[e_i, e_j]
            k+=1
    return CondensedMatrix(inner_data)
            
            