function define_base_parameters(){
    return {
                "control":{
                    "number_of_processors":undefined,
                    "algorithm_scheduler_sleep_time":undefined,
                    "scoring_scheduler_sleep_time":undefined
                },
                "workspace":{
                             "base": undefined,
                },
                
                "matrix":    {
                            "matrix_path": undefined,
                            "store_matrix_path" : undefined,
                            "pdb1": undefined,
                            "pdb2": undefined,
                            "fit_selection" : undefined,
                            "rmsd_selection": undefined
                },
                
                "algorithms":   {
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
                              "query_types":[],
                              "evaluation_criteria":[]
                },
                "results":{
                           
                }
    }
}
