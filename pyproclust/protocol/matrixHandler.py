'''
Created on 19/09/2012

@author: victor
'''
import pyproclust.tools.commonTools as common
import pyRMSD.RMSD
from pyRMSD.condensedMatrix import CondensedMatrix #@UnresolvedImport
import pickle
import numpy
import scipy.stats

class MatrixHandler(object):
    
    def __init__(self,statistics_folder=None):
        self.distance_matrix = None
        self.min_dist = None
        self.max_dist = None
        self.mean_dist = None
        self.statistics_folder = statistics_folder

    def createMatrix(self, coordsets):
        common.print_and_flush("Calculating matrix...")
        rmsd = pyRMSD.RMSD.calculateRMSDCondensedMatrix(coordsets, "OMP_CALCULATOR")
        self.distance_matrix = CondensedMatrix(rmsd)
        self.__calculate_statistics()
        self.__save_statistics()
        common.print_and_flush(" Done\n")
        
    def saveMatrix(self, matrix_file_without_extension):
        common.print_and_flush("Writing matrix data (in "+matrix_file_without_extension+".bin) ...")
        output_handler = open(matrix_file_without_extension+".bin",'w')
        pickle.dump(self.distance_matrix.get_data(),output_handler)
        output_handler.close()
        common.print_and_flush(" Done\n")
    
    def loadMatrix(self, matrix_file_without_extension):
        common.print_and_flush("Loading matrix data from "+matrix_file_without_extension+".bin ...")
        input_file_handler = open(matrix_file_without_extension+".bin","r")
        rmsd_data = pickle.load(input_file_handler)
        input_file_handler.close()
        self.distance_matrix = CondensedMatrix(list(rmsd_data))
        self.__calculate_statistics()
        self.__save_statistics()
        common.print_and_flush(" Done\n")
    
    def __calculate_statistics(self):
        data = self.distance_matrix.get_data()
        self.min_dist = numpy.min(data)
        self.max_dist = numpy.max(data)
        self.mean_dist = numpy.mean(data)
        self.std = numpy.std(data)
        self.skew = scipy.stats.skew(data)
        self.kurtosis = scipy.stats.kurtosis(data)
        
    def __save_statistics(self):
        if self.statistics_folder!=None:
            file_handler = open(self.statistics_folder+"/"+"statistics.txt","w")
            file_handler.write( "--------------\n")
            file_handler.write( "Matrix values\n") 
            file_handler.write( "-------------\n")
            file_handler.write( "Minimum:      %.5f\n"%(self.min_dist ))
            file_handler.write( "Maximum:      %.5f\n"%(self.max_dist ))
            file_handler.write( "Mean:         %.5f\n"%(self.mean_dist ))
            file_handler.write( "Std. Dev.:    %.5f\n"%(self.std))
            file_handler.write( "Skewness:     %.5f\n"%(self.skew))
            file_handler.write( "Kurtosis:     %.5f\n"%(self.kurtosis))
            file_handler.close()