function load_text_resource_with_ajax(resource){
    text_resource = ""
    
    $.ajax({
              url: resource,
              type: "GET",
              dataType: "text",
              async: false,
              
              complete: function(jqXHR, textStatus){
                  text_resource = jqXHR.responseText;
              },
              
              error:function( jqXHR, textStatus, errorThrown ){
                  alert( "Request failed: " + textStatus+". Is the server working?" );
              }
            });
            
    return text_resource;
}
