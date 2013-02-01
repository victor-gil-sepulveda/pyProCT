/*
    Function that parses webpage's controls' looking for needed parameters.
*/
function gather_params(clustering_algorithm_fields){
        
        var params = define_base_parameters();
        var param_types = define_base_parameter_types();
        
        // Fill general params
        for (var params_section in params){
            for (var params_subsection in params[params_section]){
                var param_name = params_subsection;
                var type = undefined;
                
                if($("#"+param_name).length != 0){
                    try{
                        type = param_types[params_section][params_subsection];
                    }
                    catch(error){
                        type = undefined;
                    }
                    params[params_section][params_subsection] =  get_value_of("#"+param_name, type);
                }
                else{
                    console.log(param_name+" does not exist")
                }
            }
        }

        var warnings = [] //check_regular_params(params);
        
        // Fill algorithm params
        var algorithms_section = params["algorithms"]
        for (var algorithm in algorithms_section){
                   console.log(algorithm)     
            var clustering_params = params["algorithms"][algorithm];
            var field = clustering_algorithm_fields[algorithm];
            if(!field.attr("hidden")){
                clustering_params["use"] = true;
                clustering_params["auto"] = get_value_of(field.find("#guess_params_"+algorithm));
                if(!clustering_params["auto"]){
                    try{
                        clustering_params["parameters"] = get_algorithm_parameter_parsers()[algorithm](field);
                    }
                    catch(error_message){
                        throw {message:error_message,field:field};
                    }
                }
            }
            else{
                clustering_params["use"] = false;
            }
        }
        
        console.log(params)
        
        return {"parameters":params, "warnings":warnings};
}
/**
 *
 *
 **/
function check_regular_params(params, params_description){
    var warnings = [];
    
    // /home/victor/Escritorio/mare-cuda
    
    // Check that: the root folder is not empty and exists.
    try{
        var file_check = check_path(params["workspace"]["base"], "exists", "isdir", "You must specify a project root folder.", $("#global_options_field"));
        
        if(file_check["exists"]){
             warnings.push("The folder '"+params["workspace"]["base"]+"' already exists. An execution over this\
             folder will erase its contents.");
        }
    }
    catch(error_message){
        if(error_message.existence != undefined &&  !error_message.existence){
            if(!create_folder(params["workspace"]["base"])["done"]){
                throw {message:"Impossible to create folder on that location("+params["workspace"]["base"]+").", 
                        field:$("#global_options_field")};
            }
        }
    }
    
    // Check that: if we have to load a matrix, the path is not empty and file exists.
    if(params['matrix']['creation'] == "load"){
        check_path(params['matrix']['matrix_path'], "exists", "isfile", "You must specify a location from where to load the matrix.", $("#matrix_calculation"));
    }
    
    // Check that: if we have to save a matrix, the path is not empty and file doesn't exist.
    if(params['matrix']['save_matrix']){
        check_path(params['matrix']['store_matrix_path'], "not exists", "isfile", "You must specify a location to save the matrix.", $("#matrix_calculation"));
    }
    
    // Check that: trajectories are not empty and exist
    var pdb1 = params['matrix']['pdb1'];
    var pdb2 = params['matrix']['pdb2'];
    
    check_path(pdb1, "exists", "isfile", "You must specify a trajectory file to load.", $("#datos_trayectorias"));
    console.log(params['matrix']['action'])
    if(params['matrix']['action'] == "comparison"){
        check_path(pdb2, "exists", "isfile", "When comparing trajectories you must specify a second trajectory.", $("#datos_trayectorias"));
    }
    
    // Check that: numerical values are not NaN and > 0
    var numerical = {"control": ["number_of_processors","algorithm_scheduler_sleep_time",
                                "scoring_scheduler_sleep_time"],
                      "evaluation":["maximum_noise","minimum_cluster_size",
                                "maximum_clusters","minimum_clusters"]
    };
    
    for (var params_section in numerical){
        for (var i = 0; i < numerical[params_section].length; i++){
            var numerical_param_value = params[params_section][numerical[params_section][i]]; 
            //console.log(params_section+" "+numerical_param_value)
            if(isNaN(numerical_param_value) || numerical_param_value <= 0){
                throw {message:"Field "+numerical[params_section][i]+
                                " is not a valid numerical value.", 
                                field:$(params_section)};
            }
        }
    }

    // Check that we have at least one evaluation criteria
    console.log(params['evaluation']['evaluation_criteria']);
    if($.isEmptyObject(params['evaluation']['evaluation_criteria'])){
        throw {message:"You need to define at least one criteria to select the best clustering.", 
                                field:$("#best_clustering_field")};
    }
    
    return warnings;
}
/**
 *
 *
 **/
function check_path( path, check_if, target_property, empty_message, related_field){
    if(path == ""){
        throw{
                message:empty_message, 
                field:related_field
              };
    }
    else{
        var file_check  = file_exists(path);
        
        var exists_message = "";
        if(check_if == "exists"){
            if (!file_check["exists"]){
                exists_message = "The file"+path+" does not exist.";
            }
        }
        
        if(check_if == "not exists"){
            if (file_check["exists"]){
                exists_message = "The file"+path+" already exists.";
            }
        }
        
        if(!file_check["exists"]){
            throw {message:"This file ("+path+") does not exist.", 
                    field:related_field, existence:file_check["exists"]};
        }
        
        if(!file_check[target_property]){
            var message = "This location ("+path+") is pointing to a folder."
            if(target_property == "isfile"){
                message += " It must be a file.";
            }
            if(target_property == "isdir"){
                message += " It must be a directory.";
            }
            throw {message:message, field:related_field};
        }
        
        return file_check;
    }
    
}

