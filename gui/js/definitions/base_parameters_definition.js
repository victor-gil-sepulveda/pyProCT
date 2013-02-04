function define_base_parameters(){
    return $.parseJSON(load_text_resource_with_ajax("/shared_definitions/base_parameters.json"));
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
