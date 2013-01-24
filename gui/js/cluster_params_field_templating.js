
function generate_all_fields(clustering_algorithm_fields, clustering_algorithm_ids){
    for (clustering_algorithm_name in clustering_algorithm_fields){
        clustering_algorithm_id = clustering_algorithm_ids[clustering_algorithm_name];
        clustering_algorithm_fields[clustering_algorithm_name] = {"contents":create_clustering_algorithm_field(clustering_algorithm_name,clustering_algorithm_id),
                                                                 "hidden":true};
    }
}

function create_clustering_algorithm_field(clustering_algorithm_name, clustering_algorithm_id){
    parameters_field = jQuery("<div>",{
                            class:"clustering_params_field",
                            hidden:'true',
                            style :'border:2px solid;',
                            id: clustering_algorithm_id+"_params_field"
                        });
    parameters_field.append(get_contents_for(clustering_algorithm_name, clustering_algorithm_id));
    
    
    return parameters_field;
}

function prepare_field_interactive_buttons(clustering_algorithm, clustering_algorithm_ids){
    clustering_algorithm_id = clustering_algorithm_ids[clustering_algorithm];
    
    // Prepare the checkbox
    toggable_checkbox(  "#guess_params_"+clustering_algorithm_id, 
                        "#params_field_"+clustering_algorithm_id, 
                        false);

    // Prepare the closing button
    $("#closing_button_"+clustering_algorithm_id).click(function(){
                                remove_clustering_algorithm_field(clustering_algorithm,
                                                                  clustering_algorithm_fields);
                                });
}

function add_clustering_algorithm_field(clustering_algorithm,
                                        clustering_algorithm_fields, 
                                        clustering_algorithm_ids){
    
    if(clustering_algorithm_fields[clustering_algorithm]["hidden"]){
        console.log("Adding "+clustering_algorithm);
        clustering_algorithm_fields[clustering_algorithm]["hidden"] = false;
        clustering_algorithm_fields[clustering_algorithm]["contents"].appendTo($("#"+clustering_algorithm_id+"_field"));
        clustering_algorithm_fields[clustering_algorithm]["contents"].show("slide");
        
        prepare_field_interactive_buttons(clustering_algorithm, clustering_algorithm_ids)
    }
}

function remove_clustering_algorithm_field(clustering_algorithm,
                                            clustering_algorithm_fields){
    
    if(!clustering_algorithm_fields[clustering_algorithm]["hidden"]){
        console.log("Removing "+clustering_algorithm);
        clustering_algorithm_fields[clustering_algorithm]["contents"].hide({effect:"slide",
            complete:function(){
                clustering_algorithm_fields[clustering_algorithm]["hidden"] = true;
                clustering_algorithm_fields[clustering_algorithm]["contents"].detach();
            }
        });
    }
}

function wrap_with_basic_clustering_algorithm_field(this_contents, clustering_algorithm, clustering_algorithm_id){
    return "\
            <div class ='algorithm_tittle'>"+clustering_algorithm+"<div>\
            <div class='closing_button' id='closing_button_"+clustering_algorithm_id+"'> x </div> \
            <table border='0'>\
                <tbody>\
                    <tr>\
                        <td>\
                            Let the software guess the best parameters \
                        </td>\
                        <td>\
                            <input checked='checked' id='guess_params_"+clustering_algorithm_id+"' type='checkbox'>\
                        </td>\
                    </tr>\
                </tbody>\
            </table> <div hidden = 'true' id = 'params_field_"+clustering_algorithm_id+"'>"                  
            +this_contents + "</div>"
}

function get_contents_for(clustering_algorithm, clustering_algorithm_id){
    contents = ""
    switch(clustering_algorithm){
        case "GROMOS Algorithm":
        case "Hierarchical Algorithm":
            contents ="\
                <table border='0'>\
                    <tbody>\
                        <tr>\
                            <td>\
                                Cutoff (list) : \
                            </td>\
                            <td>\
                                <input type='text'>\
                            </td>\
                        </tr>\
                    </tbody>\
                </table>";
            break;
            
        case "K-Medoids Algorithm":
            contents ="\
                <table border='0'>\
                    <tbody>\
                    <tr>\
                        <td>\
                                Seeding type : \
                            </td>\
                            <td>\
                                <select id='kmedoids_seeding_type'>\
                                    <option> RANDOM </option>\
                                    <option> EQUIDISTANT </option>\
                                    <option> GROMOS </option>\
                                </select>\
                            </td>\
                        </tr>\
                        <tr>\
                            <td>\
                                Number of clusters (list) : \
                            </td>\
                            <td>\
                                <input type='text'>\
                            </td>\
                        </tr>\
                    </tbody>\
                </table>";
            break;
            
        case "DBSCAN Algorithm":
            contents ="\
                <table border='0'>\
                    <tbody>\
                        <tr>\
                            <td>\
                                Eps (list) : \
                            </td>\
                            <td>\
                                <input type='text'>\
                            </td>\
                        </tr>\
                        <tr>\
                            <td>\
                                Minpts (list) : \
                            </td>\
                            <td>\
                                <input type='text'>\
                            </td>\
                        </tr>\
                    </tbody>\
                </table>";
            break;
        
        case "Spectral Algorithm":
            contents ="\
                <table border='0'>\
                    <tbody>\
                        <tr>\
                            <td>\
                                Sigma : \
                            </td>\
                            <td>\
                                <input type='text'>\
                            </td>\
                        </tr>\
                        <tr>\
                            <td>\
                                Number of clusters (list) : \
                            </td>\
                            <td>\
                                <input type='text'>\
                            </td>\
                        </tr>\
                    </tbody>\
                </table>";
            break;
        
        case "Random Algorithm":
            contents ="\
                <table border='0'>\
                    <tbody>\
                        <tr>\
                            <td>\
                                Number of clusters (list) : \
                            </td>\
                            <td>\
                                <input type='number'>\
                            </td>\
                        </tr>\
                    </tbody>\
                </table>";
            break;
        
        default:
            break;
        }
      
   return wrap_with_basic_clustering_algorithm_field(contents, clustering_algorithm, clustering_algorithm_id);
}
