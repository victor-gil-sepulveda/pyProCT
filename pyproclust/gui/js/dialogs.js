function create_common_dialog(type, contents, field_to_highlight, on_close_extra, parameters){
    var symbol = "";
    var title = "";
    if(type=="warning"){
        symbol = "<span class='ui-icon ui-icon ui-icon-alert' style='float: left; margin: 0 7px 50px 0;'></span>";
        title = "Warning";
    }
    
    if(type=="error"){
        symbol = "<span class='ui-icon ui-icon ui-icon-error' style='float: left; margin: 0 7px 50px 0;'></span>";
        title = "Error";
    }
    
    $("<div >", {title: title, id:'common_dialog'})
    // Add contents to the dialog
    .append(symbol+"<div class='modal_dialog'>"+
    
    contents+
    "</div>")
    // Set up dialog
    .dialog({modal:true, 
            autoResize:true,
            width:'auto',
            close: function( event, ui ){
                 on_close_extra();
                 $(this).dialog("destroy");
            },
            buttons: [{ text: "Ok",
                        click: function() { 
                            if(on_close_extra != undefined){
                                if(parameters != undefined){
                                    on_close_extra()(parameters);
                                }
                                else{
                                    on_close_extra()();
                                }
                            }
                            
                            $(this).dialog("destroy");
                            
                            if(field_to_highlight != undefined){
                                field_to_highlight.effect({effect:"highlight",duration:2000});
                            }
                        }
                             
                      }]
            });
}

