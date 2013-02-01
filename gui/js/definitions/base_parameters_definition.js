function define_base_parameters(){
    return {
                "control":{
                    "number_of_processors":undefined,
                    "algorithm_scheduler_sleep_time":undefined,
                    "scoring_scheduler_sleep_time":undefined
                },
                "workspace":{
                             "base": undefined,
                             "results": "/results",
                             "tmp": "/tmp",
                             "clusterings": "/clusterings",
                             "refinement_base" : "refinement"
                },
                
                "matrix":{
                            "creation": undefined,
                            "matrix_path": undefined,
                            "save_matrix": undefined,
                            "store_matrix_path" : undefined,
                            "action": undefined,
                            "pdb1": undefined,
                            "pdb2": undefined,
                            "fit_selection" : undefined,
                            "rmsd_selection": undefined
                },
                
                "algorithms":{
                              "gromos":{
                                       "use": undefined,
                                       "auto": undefined,
                                       "parameters":[]
                              },
                              "spectral":{
                                       "use": undefined,
                                       "auto": undefined,
                                       "sigma": [],
                                       "parameters":[]
                              },
                              "dbscan":{
                                       "use": undefined,
                                       "auto": undefined,
                                       "parameters":[]
                              },
                              "hierarchical":{
                                       "use": undefined,
                                       "auto": undefined,
                                       "parameters":[]
                              },
                              "kmedoids":{
                                       "use": undefined,
                                       "auto": undefined,
                                       "parameters":[]
                              },
                              "random":{
                                       "use": undefined,
                                       "auto": undefined,
                                       "parameters":[]
                              }
                },
                "evaluation":{
                              "maximum_noise": undefined,
                              "minimum_clusters": undefined,
                              "maximum_clusters": undefined,
                              "minimum_cluster_size":undefined,
                              "query_types":[],
                              "evaluation_criteria":[]
                },
                "refinement":{
                              "use": false,
                              "step": undefined,
                              "minimum_clusters":undefined,
                              "maximum_clusters":undefined
                },
                "results":{
                           
                }
    };
}

function define_base_parameter_types(){
    return {
                "control":{
                    "number_of_processors":"int",
                    "algorithm_scheduler_sleep_time":"int",
                    "scoring_scheduler_sleep_time":"int"
                },
                
                "matrix":{
                         "creation": "radio",
                         "save_matrix": "checkbox",
                         "action": "radio",
                         "pdb1": "text",
                         "pdb2": "text",
                         "matrix_path": "text",
                         "store_matrix_path" : "text",
                },
                "evaluation":{
                              "maximum_noise": "int",
                              "minimum_clusters": "int",
                              "maximum_clusters": "int",
                              "minimum_cluster_size":"int",
                              "query_types":"tags",
                              "evaluation_criteria":"tags:criteria"
                }
    };
}
