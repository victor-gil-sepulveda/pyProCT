/**
 *
 **/
function set_up_run(clustering_algorithm_fields){
    return function(){
        // A 'gathering_result' contains a description of warnings and parameters.
        // If it fails, the operation is aborted and an error dialog is shown.
        var gathering_result = undefined;
        try{
            var gathering_result = gather_params(clustering_algorithm_fields);
        }
        catch(exception){
            create_common_dialog("error", exception.message, exception.field);
            // End the game here.
            return;
        }
        
        // If there was not a fatal error, we will show the informative warnings if any.
        var warnings = gathering_result["warnings"];
        var final_warning_text = "";
        for(var i =0; i < warnings.length;i++){
            final_warning_text += "<br>&#9;-"+warnings[i];
        }
        console.log(final_warning_text)
        // At this point we can generate the script
        var parameters = gathering_result["parameters"];
        if(final_warning_text != ""){
            create_common_dialog("warning", 
            "The following warnings were generated:"+
            final_warning_text+
            "<br>The script will be generated, but its execution can give problems.", undefined, ajax_save_script, parameters);
        }
        else{
            ajax_run_script()(parameters);
        }
    }
}

/**
 *
 **/
function ajax_run_script(){
        return function(parameters){
        $.ajax({
              url: "/run",
              type: "POST",
              data: JSON.stringify(parameters),
              dataType: "text",
              
              complete: function(jqXHR, textStatus){
                  console.log(jqXHR.responseText)
                  var my_response_object =  $.parseJSON(jqXHR.responseText);
                  start_monitoring_run();
              },
              
              error:function( jqXHR, textStatus, errorThrown ){
                  alert( "Request failed: " + textStatus+". Is the server working?" );
              }
            });
        }
}       


function start_monitoring_run(){

}
