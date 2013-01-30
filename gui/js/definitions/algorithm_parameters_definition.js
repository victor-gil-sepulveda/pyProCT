function get_algorithm_parameters_definition(){
    return {
            "gromos":[{
                        label:"Cutoff (list) :",
                        control_id: "cutoff_list",
                        control_type: "text",
                      }
            ],
            "hierarchical":[
                            {
                                label:"Cutoff (list) :",
                                control_id: "cutoff_list",
                                control_type: "text"
                            },
                            {
                                label:"Method :",
                                control_id: "method_type",
                                control_type: "select",
                                options: ["COMPLETE","SINGLE","AVERAGE"]
                            }
            ],
            "kmedoids":[
                        {
                            label:"Seeding type :",
                            control_id: "kmedoids_seeding_type",
                            control_type: "select",
                            options: ["RANDOM", "EQUIDISTANT", "GROMOS" ]
                        },
                        {
                            label:"Number of clusters (list) : ",
                            control_id: "number_of_clusters",
                            control_type: "text"
                        }
            ],
            "dbscan":[
                        {
                            label:"Eps (list) : ",
                            control_id: "eps_list",
                            control_type: "text"
                        },
                        {
                            label:"Minpts (list) : ",
                            control_id: "minpts_list",
                            control_type: "text"
                        }
            ],
            "spectral":[
                        {
                            label:"Sigma (list) : ",
                            control_id: "sigma",
                            control_type: "text"
                        },
                        {
                            label:"Number of clusters (list) : ",
                            control_id: "number_of_clusters",
                            control_type: "text"
                        }
            ],
            "random":[
                      {
                        label:"Number of clusters (list) : ",
                        control_id: "number_of_clusters",
                        control_type: "text"
                      }
            ]
        };
}
