'''
Created on 20/02/2012

@author: victor
'''

import json
import os.path
import validation.bidimensional.datasets as data
import validation.bidimensional.validationTools as vtools
from pyproct.tools.scriptTools import create_directory
from pyproct.driver.handlers.matrix.matrixHandler import MatrixHandler
from pyproct.driver.parameters import ProtocolParameters
from pyproct.driver.observer.observer import Observer
from pyproct.driver.driver import Driver
from pyproct.tools.commonTools import convert_to_utf8
from pyproct.clustering.clustering import Clustering

if __name__ == '__main__':
    create_directory("./clustering_images")
    create_directory("./matrices")
    create_directory("./tmp")
    condensed_matrices, all_observations = vtools.create_matrices(data)
    # Saving matrices
    for dataset_name in data.all_datasets:
        handler = MatrixHandler({"method":"load"})
        handler.distance_matrix = condensed_matrices[dataset_name]
        handler.save_matrix("./matrices/%s"%dataset_name)

    # Run pyProCT for each of them
    base_script = "".join(open("base_script.json","r").readlines())
    for dataset_name in data.all_datasets: #["spaeth_06"]:#
        print dataset_name
        # Change placeholders
        script_str = base_script%(os.path.abspath("./tmp/%s"%dataset_name),"./matrices/%s"%dataset_name)
        parameters = ProtocolParameters.get_params_from_json(script_str)
        # And change another hypothesis stuff
        parameters["clustering"]["evaluation"]["maximum_noise"] = data.noise[dataset_name]
        parameters["clustering"]["evaluation"]["minimum_cluster_size"] = data.minsize[dataset_name]
        parameters["clustering"]["evaluation"]["minimum_clusters"] = data.num_cluster_ranges[dataset_name][0]
        parameters["clustering"]["evaluation"]["maximum_clusters"] = data.num_cluster_ranges[dataset_name][1]
        print parameters["clustering"]["evaluation"]["minimum_clusters"], parameters["clustering"]["evaluation"]["maximum_clusters"]
        if dataset_name in data.criteria:
            parameters["clustering"]["evaluation"]["evaluation_criteria"] = data.criteria[dataset_name]
        else:
            parameters["clustering"]["evaluation"]["evaluation_criteria"] = data.criteria["default"]
        Driver(Observer()).run(parameters)

    for dataset_name in data.all_datasets:
        results_file = os.path.join(os.path.abspath("./tmp/%s"%dataset_name),"results/results.json")
        results = convert_to_utf8(json.loads(open(results_file).read()))
        best = results["best_clustering"]
        clustering = Clustering.from_dic(results["selected"][best]["clustering"])
        vtools.show_2D_dataset_clusters(all_observations[dataset_name],
                                        clustering,
                                        scale = 20,
                                        margin = 20).save("clustering_images/%s.jpg"%dataset_name,
                                                 "JPEG")
        print dataset_name,results["selected"][best]["type"],results["selected"][best]["clustering"]["number_of_clusters"], results["selected"][best]["evaluation"]["Noise level"],#results["selected"][best]["parameters"]
        # look for the best criteria
        criteria_scores = []
        for criteria in results["scores"]:
            criteria_scores.append((results["scores"][criteria][best],criteria))
        print criteria_scores

    print "\nDone"