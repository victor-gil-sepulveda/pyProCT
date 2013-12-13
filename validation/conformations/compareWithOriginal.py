'''
Created on 11/12/2013

@author: victor
'''
from pyproct.driver.parameters import ProtocolParameters
from pyproct.driver.driver import Driver
from pyproct.driver.observer.observer import Observer
import sys
from pyproct.driver.handlers.matrix.matrixHandler import MatrixHandler
from pyproct.tools.matrixTools import get_submatrix
from pyproct.tools.plotTools import matrixToImage

if __name__ == '__main__':
    base_script = "".join(open("base_script.json","r").readlines())
    parameters = ProtocolParameters.get_params_from_json(base_script)
    parameters["workspace"]["base"] = sys.argv[3]
    parameters["global"]["pdbs"] = [sys.argv[1], sys.argv[2]]
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

        matrixToImage(get_submatrix(handler.create_matrix(None), range(9,18)), parameters["workspace"]["base"] +"/submatrix.png")