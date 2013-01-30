/*
    Returns the parameters parser for a given algorithm.
*/
function get_algorithm_parameter_parsers(){
    return {   
                "gromos":parse_gromos_parameters, 
                "kmedoids": parse_kmedoids_parameters,
                "hierarchical":parse_hierarchical_parameters,
                "spectral":parse_spectral_parameters,
                "dbscan":parse_dbscan_parameters,
                "random":parse_random_parameters
            };
}

/*
    Creates a copy of one parameter.
*/
function clone_param(this_param){
    var cloned = {};
    for (key in this_param){
        cloned[key] = this_param[key];
    }
    return cloned;
}

/**
 *   Extends a parameter list with an attribute that can have diverse values.
 *   For instance, extending [{"a":3,"b":4},{"d":7}] with "c" = [5,6] would produce:
 *   [{"a":3,"b":4,"c":5} , {"a":3,"b":4,"c":6}, {"d":7,"c":5}, {"d":7,"c":6}]
 *
 *      
 **/
function combine_parameters(this_parameters, with_this_parameter, with_this_values){
    var extended_parameters = [];
    for(var i = 0; i < with_this_values.length; i++){
        for(var j = 0; j < this_parameters.length; j++){
           var new_param = clone_param(this_parameters[j]);
           new_param[with_this_parameter] =  with_this_values[i];
           extended_parameters.push(new_param);
        }
    }
    return extended_parameters;
}

/**
 *   Generates the parameters for one or more list of values.
 *   Example:
 *   bind_parameters({"a":[1,2],"b":[3,4]}) == [{"a":1,"b":3},{"a":2,"b":4}]
 *   bind_parameters({"a":[1,2],"b":[3,4]},{"d":2}) ==  [{"a":1,"b":3,"d":2},{"a":2,"b":4,"d":2}]
 *
 *   @param {As. Array} lists_dic is a dictionary where the key will be the name of the parameter,
 *   and the value for this key will be the list of values for this parameter.
 *
 *   @param {As. Array} singular a dictionary with values that remain constant for all parameters.
 *
 *   @returns {list} A list of parameters.      
 **/
function bind_parameters(lists_dic,singular){
    var parameters = [];
    // length
    var len = 0;
    for (var key in lists_dic) {
        len = lists_dic[key].length;
    }
    
    // Create parameters from lists
    for(var i = 0; i < len; i++){
        var parameter = {};
        for(name in lists_dic){
            parameter[name] = lists_dic[name][i]; 
        }
        parameters.push(parameter);
    }
    
    // Singular values
    if(singular != undefined){
        for(var i = 0; i < len; i++){
            for(sing in singular){
                parameters[i][sing] = singular[sing];
            }
        }
    }
    
    return parameters;
}

function parse_gromos_parameters(my_field){
    var cutoff_list = parse_list(my_field.find("#cutoff_list").val(), parseFloat);
    
    check([{ 
            "rule": list_not_empty,
            "params":{
                    "cutoff list":cutoff_list
            }  
           }
    ], "gromos");
    
    return bind_parameters({
                            "cutoff":cutoff_list
                            });
}

function parse_hierarchical_parameters(my_field){
    var cutoff_list = parse_list(my_field.find("#cutoff_list").val(), parseFloat);
    
    check([
            { "rule": list_not_empty,
              "params":{
                    "cutoff list":cutoff_list
              }  
            }
    ], "hierarchical");
    
    var method = get_value_of(my_field.find("#method_type"));
    
    return bind_parameters({
                            "cutoff":cutoff_list
                            },
                            {
                            "method":method
                            });
}

function parse_dbscan_parameters(my_field){
    var eps_list = parse_list(my_field.find("#eps_list").val(), parseFloat);
    var minpts_list = parse_list(my_field.find("#minpts_list").val(), parseFloat);
    
    check([
            { "rule": list_not_empty,
              "params":{
                    "eps list":eps_list
              }  
            },
            { "rule": list_not_empty,
              "params":{
                    "minpts list":minpts_list
              }  
            },
            { "rule": lists_have_same_length,
              "params":{
                    "eps list":eps_list,
                    "minpts list":minpts_list
              }  
            }
    ], "dbscan");
    
    return bind_parameters({
                                "eps":eps_list,
                                "minpts": minpts_list        
                            });
}

function parse_kmedoids_parameters(my_field){
    var clustering_size_list = parse_list(my_field.find("#number_of_clusters").val());
    
    check([
            { "rule": list_not_empty,
              "params":{
                    "clustering size list":clustering_size_list
              }  
            }
    ], "kmedoids");
    
    /*if(seeding_type == "gromos"){
        param["seeding_max_cutoff"] = parseInt(get_value_of(my_field.find("#kmedoids_seeding_type")));
    }*/
    
    return bind_parameters({
                             "k": clustering_size_list      
                            });
}

function parse_spectral_parameters(my_field){
    var sigma_list = parse_list(my_field.find("#sigma").val(), parseFloat);
    var k_list = parse_list(my_field.find("#number_of_clusters").val());
    
    check([
            { "rule": list_not_empty,
              "params":{
                    "sigma list":sigma_list
              }  
            }
    ], "spectral");
    
    check([
            { "rule": list_not_empty,
              "params":{
                    "k list":k_list
              }  
            }
    ], "spectral");
    
    return combine_parameters(bind_parameters({
                                                //"use_k_medoids": ,
                                                "k": k_list        
                                                }), 
                                                "sigma", sigma_list );
}

function parse_random_parameters(my_field){
    var number_of_clusters_list = parse_list(my_field.find("#number_of_clusters").val());
    
    check([{ 
            "rule": list_not_empty,
            "params":{
                    "num_clusters": number_of_clusters_list
            }  
           }
    ], "random");
    
    return bind_parameters({
                             "num_clusters": number_of_clusters_list      
                            });
}

/*
    Parses the contents of a text control holding a list of numbers description. This list can 
    have two forms:
    - Comma separated list of numbers Ex. "1, 2, 3, 4"
    - Range with this form : start, end : step  Ex. "4, 14 :2"  = "4, 6, 8, 10, 12"
*/
function parse_list(list_string, conversor){
    
    // Default value for conversor
    if (conversor == undefined){
        conversor = parseInt
    }
    
    try{
        // Remove non dot, colon digit or character 
        var list_string = list_string.replace(/[^\d,.:]+/g, '');
        
        var sequence = [];
        
        // Analyze the string    
        var parts = list_string.split(":")
        if (parts.length == 2){
            // Is the description of a range i,j :step, those can only be integers
            var range_parts = parts[0].split(",");
            if( range_parts.length != 2){
                return undefined;
            }
            var i = parseInt(range_parts[0]);
            var j = parseInt(range_parts[1]);
            if( isNaN(i) || isNaN(j)){
                return undefined;
            }
            var step = parseInt(parts[1]);
            
            for(var k = i; k<j; k+=step){
                sequence.push(k);
            }
        }
        else {
            // the list is a comma-separated sequence of numbers
            var string_sequence = parts[0].split(",");
            for(var i = 0; i < string_sequence.length; i++){
            console.log("*"+string_sequence[i]+"*")
                var value = conversor(string_sequence[i]);
                if(!isNaN(value)){
                    sequence.push(value);
                }
            }
        }
    }
    catch(error_message){
        throw "There was an error while parsing this list:["+list_string+"]. Error was: "+error_mesage
    }
    
    return sequence;
}

/*
    Helper function to get the value of an undetermined control.
*/
function get_value_of(of_this){
    switch($(of_this).attr('type')){
        case "text":
            return $(of_this).val();
        case "checkbox":
            return $(of_this).is(":checked");
        case "number":
            return $(of_this).val();
        default:
            return $(of_this).val();
    }
}

