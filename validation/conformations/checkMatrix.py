"""
Created on 29/01/2014

@author: victor
"""
import sys
from pyproct.driver.handlers.matrix.matrixHandler import MatrixHandler
import numpy

if __name__ == '__main__':
    matrix_path = sys.argv[1]
    handler = MatrixHandler({
            "method": "load",
            "parameters": {
                "path": matrix_path
            }
        })
    matrix = handler.create_matrix(None)
    print matrix.get_number_of_rows()
    data =  list(matrix.get_data())*2+[0.]*9
    print "Avg. %0.4f"%numpy.mean(data)
    print "Std. dev. %0.4f"%(numpy.std(data))

