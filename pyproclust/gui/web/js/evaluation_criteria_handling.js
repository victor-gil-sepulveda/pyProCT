/*
    Prepares and shows the evaluation criteria dialog.
*/
function criteria_creation_show_dialog(criteria_types, tag_widget_id, template){
    return function(){
        $("<div >", {title: "New Criteria",id:'criteria_creation_dialog'})
        // Add contents to the dialog
        .append(get_eval_dialog_contents(criteria_types, template))
        // Set up dialog
        .dialog({modal:true, 
                autoResize:true,
                width:'auto',
                create: function( event, ui ) {
                   $(".dialog_spinner").spinner({places:2,step:0.05});
                   $(".dialog_spinner").css({width:"35px"});
                },
                close: function( event, ui ){
                     $(this).dialog("destroy"); 
                },
                buttons: [{ text: "Discard",
                            click: function() { 
                                $(this).dialog("destroy"); 
                                } 
                          },
                          { text: "Ok",
                            click: function() { 
                                var criteria = criteria_to_string('criteria_creation_dialog', criteria_types);
                                $("#"+tag_widget_id).tagit("createTag",criteria);
                                $(this).dialog("destroy");
                                } 
                          }]
                })
    };
}
/*
    Creates the contents of the dialog (using handlebars))
*/
function get_eval_dialog_contents(criteria_list, template){
    // Gather data
    var data = {criteria:[]};
    for (var i = 0;i < criteria_list.length; i++){
        var criteria_name = criteria_list[i];
        data.criteria.push({name:criteria_name,  initial_value:0});
    }
    
    // Render it
    //var source   = $("#dialog_contents_template").html();
    var template = Handlebars.compile(template);
    return template(data);
}

/*
    Creates a string from the contents of the dialog that represents one evaluation criteria.
*/
function criteria_to_string(dialog_to_extract_data, criteria_list){
    var string_criteria = "";
    for (var i =0; i<criteria_list.length; i++){
        var criteria_name = criteria_list[i];
        // Get the spinners value, if different than 0, proceed
        var weight =  $("#"+criteria_name+"_spinner").val();
        if (weight != 0){
            // Maximize or minimize?
            var min_max = $("#"+criteria_name+"_listbox").val();
            string_criteria += min_max + " " + criteria_name + " (weigth: "+ weight + ") and"   
        }
    }
    // remove last and return
    return string_criteria.substring(0,string_criteria.length-4);
}
