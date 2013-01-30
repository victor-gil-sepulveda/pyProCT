
function set_up_param_saving(clustering_algorithm_fields){
    return function(){
                var request = $.ajax({
                  url: "save_params",
                  type: "POST",
                  data: JSON.stringify(gather_params(clustering_algorithm_fields)),
                  dataType: "text",
                  complete: function(jqXHR, textStatus) {
                      console.log(jqXHR.responseText)
                      var myo =  $.parseJSON(jqXHR.responseText);
                      alert(myo.caca)
                  },
                  error:function( jqXHR, textStatus, errorThrown ){
                      alert( "Request failed: " + textStatus );
                  }
                });
    }
}                 


/*
    Function that parses webpage's controls' looking for needed parameters.
*/
function gather_params(clustering_algorithm_fields){
        
        params = define_base_parameters();
        
        // Fill general params
        for (var params_section in params){
            console.log("*"+params_section)
            for (var params_subsection in params[params_section]){
                var param_name = params_subsection;
                
                if($("#"+param_name).length != 0){
                    params[params_section][params_subsection] =  get_value_of("#"+param_name);
                }
                else{
                    console.log(param_name+" does not exist")
                }
            }
        }
        
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
                        alert("Error: "+error_message);
                    }
                }
            }
            else{
                clustering_params["use"] = false;
            }
        }
        
        console.log(params)
        
        return params;
}
