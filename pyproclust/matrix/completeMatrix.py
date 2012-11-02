'''
Created on 12/03/2012

@author: victor
'''
import math
from pyproclust.matrix.matrixCommon import generate_row

def load_complete_matrix(file_handler):
    """ 
    Loads a matrix from a file where each line is a row.
    """
    m = CompleteDistanceMatrix()
    m_contents = m.get_data()
    for l in file_handler:
        row = []
        generate_row(l,row) 
        m_contents.append(row)
    return m


class CompleteDistanceMatrix(object):
    """
    This class is a representation of a symmetric square matrix used to store
    a distance matrix.
    """
    def __init__(self,contents = []):
        """
        Constructor.
        """
        self.contents = contents
    
    def get_row_dimension(self):
        """
        Returns the length of a row (and thus of the columns).
        """
        return len(self.contents)
    
    def get_data(self):
        """
        Returns the internal representation.
        """
        return self.contents
    
    def reflect(self):
        """
        Performs a horizontal reflection to a matrix.
        """    
        jmax = len(self.contents)-1
        #print "max", len(matrix)
        for j in range((jmax+1)/2):
            # swap
            self.contents[j], self.contents[jmax-j] = self.contents[jmax-j], self.contents[j]
    
    def save(self,file_handler):
        """
        Writes the matrix into a file.
        """
        for row in self.contents:
            for number in row:
                file_handler.write(str(number)+" ")
            file_handler.write("\n")
    
    def symmetry_error(self):
        """
        Calculates the mean error of the symmetry (so, the mean of the difference of
        expected values). 
        """
        accum_error = 0.
        erroneous_elements = 0     
        for i in range(len(self.contents)):
            for j in range(i+1,len(self.contents[0])):
                if self.contents[j][i] != self.contents[i][j]:
                    erroneous_elements = erroneous_elements +1 
                    accum_error = accum_error + math.fabs(self.contents[i][j]-self.contents[j][i])
        
        if erroneous_elements > 0:
            mean_error = accum_error / erroneous_elements
            return  mean_error
        else:
            return 0.0
    
    def validate_dimensions(self, ncols, nrows, verbose = True):
        """
        Checks the dimension of a matrix
        """
        if nrows != len(self.contents):
            if verbose:
                print "Wrong column dimension. Matrix lenght:",len(self.contents), "instead of", nrows
            return False
        else:
            i = 0
            for row in self.contents:
                if ncols != len(row):
                    if verbose:
                        print "Wrong row dimension. Row:", i,"Lenght:",len(row), "instead of", ncols
                    return False
                i = i+1
        return True
    
    def get_minimum_and_maximum(self):
        """
        Returns a tuple with the minimum and maximum value of the matrix.
        """
        max_number = 0.
        min_number = 100000.
        for row in self.contents:
            max_number = max(max_number,max(row))
            min_number = min(min_number,min(row))
        return (min_number,max_number)
    
    def normalize(self, min_max):
        """
        Normalizes the values of the matrix in the range [0,1] given a minimum and 
        maximum value of the elements in the matrix.
        """
        n_range = float(min_max[1]-min_max[0])
        for row in self.contents:
            for i in range(len(row)):
                row[i] = (row[i] - min_max[0]) / n_range 
                
    def __eq__(self,other):
        """
        Equality operator.
        """
        if self.get_row_dimension() != other.get_row_dimension():
            return False
        else:
            dim = self.get_row_dimension()
            for i in range(dim):
                for j in range(dim):
                    if(self.get_data()[i][j] != other.get_data()[i][j]):
                        return False 
        return True

