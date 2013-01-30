/*
    Pregenerates all the clustering algorithm parameters fields. we are going to hide/unhide them
    depending if the algorithm is used or not.
*/
function generate_all_fields(clustering_algorithm_titles, algorithm_parameters){
    var fields = {};
    
    for (clustering_algorithm in clustering_algorithm_titles){
        var title = clustering_algorithm_titles[clustering_algorithm];
        fields[clustering_algorithm] = create_clustering_algorithm_field(title, clustering_algorithm, algorithm_parameters);
    }
    
    return fields;
}
/*
    Creates the algorithm parameters field of the 'clustering_algorithm' algorithm.
*/
function create_clustering_algorithm_field(title, clustering_algorithm, algorithm_parameters){
    parameters_field = $("<div>",{
                            class:"clustering_params_field",
                            hidden:'true',
                            style :'border:2px solid;',
                            id: clustering_algorithm+"_params_field"
                        });
    
    parameters_field.html(get_contents_for(title, clustering_algorithm, algorithm_parameters));
    
    return parameters_field;
}

/*
    Uses the template to create one of the clustering algorithm fields.
*/
function get_contents_for(title, clustering_algorithm, algorithm_parameters){
    
    var data = {
        "id": clustering_algorithm,
        "title": title,
        "properties":algorithm_parameters[clustering_algorithm]
    }
           
    var source   = $("#basic_algo_field_template").html();
    var template = Handlebars.compile(source);
    return template(data);
}

/*
    Handler for field opening (linked to the 'add clustering' button).
*/
function add_clustering_algorithm_field(clustering_algorithm, clustering_algorithm_fields){
    
    if(clustering_algorithm_fields[clustering_algorithm].attr("hidden")){
        console.log("Adding "+clustering_algorithm);
        clustering_algorithm_fields[clustering_algorithm].appendTo($("#"+clustering_algorithm+"_field"));
        clustering_algorithm_fields[clustering_algorithm].show({
                   
                   effect:"slide",
                   
                   complete:function(){
                        clustering_algorithm_fields[clustering_algorithm].attr("hidden", false);
                        prepare_field_interactive_buttons(clustering_algorithm, clustering_algorithm_fields);
                   }
        });
        
    }
}

/*
    Helper function to add properties to the controls inside the field.
*/
function prepare_field_interactive_buttons(clustering_algorithm, clustering_algorithm_fields){
                                            
    // Prepare the checkbox
    toggable_checkbox(  "#guess_params_"+clustering_algorithm, 
                        "#params_field_"+clustering_algorithm);

    // Prepare the closing button
    $("#closing_button_"+clustering_algorithm).click(function(){
                                var field = clustering_algorithm_fields[clustering_algorithm];
                                remove_clustering_algorithm_field(field);
    });
}

/*
    Handler for field closing.
*/
function remove_clustering_algorithm_field(clustering_algorithm_field){
    
    if(!clustering_algorithm_field.attr("hidden")){
        console.log("Removing "+clustering_algorithm_field.attr("id"));
        clustering_algorithm_field.hide({
            effect:"slide",
            complete:function(){
                clustering_algorithm_field.attr("hidden",true);
                clustering_algorithm_field.remove();
            }
        });
    }
}


