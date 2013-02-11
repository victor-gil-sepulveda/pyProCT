'''
Created on 12/03/2012

@author: victor
'''
import sys
import random

def cluster_from_tuple(mytuple):
    """
    Creates a cluster from a tupple formed by a first element being the prototype and
    a second element being the list of all elements.
    """
    prototype = mytuple[0]
    all_elements = list(mytuple[1])
    all_elements.extend([prototype])
    return Cluster(prototype, all_elements)

def get_cluster_sizes(clusters):
    """
    Calculates all the sizes of a clusters list and returns it in a tuple in which
    the first element is the total number of elements, and the second a list containing
    the size of each cluster (maintaining the ordering).
    """
    total_elements = 0
    cluster_sizes = []
    for c in clusters:
        size = c.get_size()
        total_elements = total_elements + size
        cluster_sizes.append(size)
    return total_elements,cluster_sizes

def gen_clusters_from_class_list(group_list,skip_list=[]):
    """
    Generates the clusters that describe a group list. A group list for a N elements clustering
    is defined as the list of N elements with the number of cluster to which each cluster
    belongs. Example: for 4 elements [1,2,3,4] a possible group list would be: [2,1,2,1] which
    means that element 0 and 2 belong to cluster 2 and the others to cluster 2. As it's not possible
    to define a centroid or medioid. ATENTION: the first element of the cluster will be defined as the
    centroid/medoid.
    """
    dic_clusters = {}
    for i in range(len(group_list)):
        if not group_list[i] in skip_list:
            if group_list[i] in dic_clusters.keys():
                dic_clusters[group_list[i]].append(i)
            else:
                dic_clusters[group_list[i]] = [i]
    clusters = []
    for k in dic_clusters.keys():
        clusters.append(Cluster(dic_clusters[k][0],dic_clusters[k]))
    return clusters

class Cluster(object):
    '''
    A cluster object is defined a group of elements which have one or more characteristics in common
    and one element which is the most representative element of the cluster. 
    '''
    most_representative_element = None
    all_elements = []
    distance_calculator = []
    
    def __init__(self, prototype , elements):
        '''
        TODO: DOCUMENT
        '''
        self.set_elements(elements)
        try:
            self.set_prototype(prototype)
        except TypeError:
            raise

    def set_prototype(self,this_one):
        '''
        Adds a representative element which must already be inside the
        internal elements list. 
        '''
        if this_one == None:
            self.prototype = None
        else:
            if this_one in self.all_elements:
                self.prototype = this_one
            else:
                raise TypeError("[Error in Cluster::set_prototype] the prototype is not in the elements list.")
    
    def set_elements(self,elements):
        self.all_elements = elements
    
    def get_size(self):
        """
        Returns the size of the cluster (which is indeed the size of its elements list)
        """
        return len(self.all_elements)

    def __eq__(self, other):
        """
        Checks whether two clusters are equal or not. Returns True or False depending
        on it :P
        """
        if(self.get_size() != other.get_size()):
            return False
        else:
            elements = sorted(self.all_elements)
            other_elements = sorted(other.all_elements)
            for i in range(len(elements)):
                if elements[i] != other_elements[i]:
                    return False
            return True

    def __str__(self):
        return "["+str(self.prototype)+str(self.all_elements)+"]"
    
    def __getitem__(self, index):
        return self.all_elements[index]
    
    def calculate_biased_medoid(self,condensed_distance_matrix,elements_into_account):
        """
        Calculates the medoid (element with minimal distance to all other objects) of the 
        elements of the cluster which are in elements_into_account. 
        """
        all_elems_set = set(self.all_elements)
        accountable_set = set(elements_into_account)
        
        # Check that elements_into_account is a subset of all_elements
        elem_inters = all_elems_set.intersection(accountable_set)
        if len(elem_inters) != len(elements_into_account):
            print "[ERROR Cluster::calculate_biased_medoid] 'elements_into_account' is not a subset of the"  
            print "elements of this cluster."
            exit()
        
        if len(elements_into_account) == 0:
            return None
        
        min_dist_pair = (sys.maxint, elements_into_account[0])
        for ei in elements_into_account:
            # Calculate distances for this vs all the others
            # Note that for comparing, the mean is not required,as
            # all have the same amount of elements
            summed_distance = 0
            for ej in elements_into_account:
                summed_distance = summed_distance +condensed_distance_matrix[ei,ej]
            min_dist_pair = min(min_dist_pair,(summed_distance,ei))
        
        medoid_element = min_dist_pair[1]
        
        return medoid_element
 
    def calculate_medoid(self,condensed_distance_matrix):
        """
        Calculates the medoid for all_elements of the cluster.
        """
        return self.calculate_biased_medoid(condensed_distance_matrix,self.all_elements)
    
    def get_random_sample(self, n, rand_seed = None):
        """
        Returns a random sample of the elements.
        
        @param n: Number of random elements to get.
        
        @param rand_seed: Seed for the random package. Used for testing (repeteability)
        
        @return: A random sample of the cluster elements.
        """
        if not rand_seed is None:
            random.seed(rand_seed)
        temporary_list = list(self.all_elements)
        random.shuffle(temporary_list)
        return temporary_list[0:n]
        