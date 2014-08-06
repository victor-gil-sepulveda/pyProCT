"""
Created on 11/12/2013

@author: victor
"""
from pyproct.driver.parameters import ProtocolParameters
from pyproct.driver.driver import Driver
from pyproct.driver.observer.observer import Observer
import sys
from pyproct.driver.handlers.matrix.matrixHandler import MatrixHandler
from pyproct.tools.matrixTools import get_submatrix
from pyproct.tools.plotTools import matrixToImage
from itertools import product
import numpy
from pyproct.tools.pdbTools import get_number_of_frames

if __name__ == '__main__':
    base_script = "".join(open("base_script.json","r").readlines())
    parameters = ProtocolParameters.get_params_from_json(base_script)
    parameters["global"]["workspace"]["base"] = sys.argv[3]
    parameters["data"]["files"] = [sys.argv[1], sys.argv[2]]

    frames_ini = get_number_of_frames(sys.argv[1])
    frames_proto = get_number_of_frames(sys.argv[2])
    print sys.argv[1],"->",frames_ini
    print sys.argv[2],"->",frames_proto

    try:
        Driver(Observer()).run(parameters)
    except SystemExit:
        # Expected improductive search
        # Load again the matrix
        handler = MatrixHandler({
            "method": "load",
            "parameters": {
                "path": parameters["workspace"]["base"]+"/matrix/matrix"
            }
        })
        matrix = handler.create_matrix(None)
        submatrix = get_submatrix(matrix, range(frames_ini,frames_ini+frames_proto))
        matrixToImage(submatrix, parameters["workspace"]["base"] +"/submatrix.png")
        print "Original mean:",get_submatrix(matrix, range(0,frames_ini)).calculateMean()
        values = []
        for i in range(0,frames_ini):
            for j in range(frames_ini,frames_ini+frames_proto):
                values.append((handler.distance_matrix[i,j],i,j-frames_ini))
        for d,i,j in sorted(values):
            print "%d %d %.2f"% (i,j,d)

        print "Combined mean:", numpy.mean(values)
