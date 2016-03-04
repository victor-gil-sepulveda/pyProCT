"""
Created on Mar 4, 2016

@author: victor
"""

class PairwiseRMSDProvider(object):
    
    def __init__(self,  data_handler):
        self.data_handler = data_handler
        
    def __getitem__(self, key):
        pass
    
    def __setitem__(self, key, value):
        pass
    
    def __len__(self):
  
        
#         "get_number_of_rows", (PyCFunction) condensedMatrix_get_number_of_rows, METH_NOARGS,PyDoc_STR("description")},
#     {"get_data", (PyCFunction) condensedMatrix_get_data, METH_NOARGS,PyDoc_STR("description")},
# 
#     // Statistics
#     {"recalculateStatistics", (PyCFunction) condensedMatrix_calculate_statistics, METH_NOARGS,PyDoc_STR("description")},
#     {"calculateMean",         (PyCFunction) condensedMatrix_get_mean, METH_NOARGS,PyDoc_STR("description")},
#     {"calculateVariance",     (PyCFunction) condensedMatrix_get_variance, METH_NOARGS,PyDoc_STR("description")},
#     {"calculateSkewness",     (PyCFunction) condensedMatrix_get_skewness, METH_NOARGS,PyDoc_STR("description")},
#     {"calculateKurtosis",     (PyCFunction) condensedMatrix_get_kurtosis, METH_NOARGS,PyDoc_STR("description")},
#     {"calculateMax",         (PyCFunction) condensedMatrix_get_max, METH_NOARGS,PyDoc_STR("description")},
#     {"calculateMin",         (PyCFunction) condensedMatrix_get_min, METH_NOARGS,PyDoc_STR("description")},
# 
#     // Matrix as graph
#     {"get_neighbors_for_node", (PyCFunction)condensedMatrix_get_neighbors_for_node, METH_VARARGS,PyDoc_STR("description")},
#     {"choose_node_with_higher_cardinality", (PyCFunction)condensedMatrix_choose_node_with_higher_cardinality, METH_VARARGS,PyDoc_STR("description")},
#     {"element_neighbors_within_radius"pass

# Py_ssize_t condensedMatrix_length(CondensedMatrix *self){
#     return self->row_length;
# }