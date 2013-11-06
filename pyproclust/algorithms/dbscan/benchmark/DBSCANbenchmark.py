import time
import random
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproclust.algorithms.dbscan.dbscanTools import k_dist

print "Creating data..."
t0 = time.time()
row_size = 5000
matrix_elem_size = row_size*(row_size-1)/2
contents = numpy.array(random.sample(xrange(matrix_elem_size+1),matrix_elem_size))
contents_range = contents.max()-contents.min()
contents = (contents-contents.mean())/contents_range
matrix = CondensedMatrix(contents)
print "It took",time.time() - t0, "seconds to create the matrix."

times = []
for i in range (5):
    t0 = time.time()
    k_list = range(5)
    k_dist(k_list, matrix)
    times.append(time.time() - t0)
times = numpy.array(times)
print "It took %.3f (%.3f) seconds to calculate kdists."%(times.mean(),times.std())

# 7.594 (0.029) -> base version
# 7.667 (0.057) -> without reshaping (using transpose)
# 7.326 (0.067) -> changing mergesort to quicksort
# 7.288 (0.021) -> with both quicksorts