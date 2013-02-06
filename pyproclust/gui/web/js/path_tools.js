/**
 *
 *
 *
 *
 **/
function file_exists(location){
    var response = undefined;
    $.ajax({
      url: "/file_exists",
      type: "POST",
      data: JSON.stringify({'location':location}),
      dataType: "text",
      async: false,
      
      complete: function(jqXHR, textStatus) {
          console.log(jqXHR.responseText)
          response =  $.parseJSON(jqXHR.responseText);
      },
      
      error:function( jqXHR, textStatus, errorThrown ){
          alert( "Request failed: " + textStatus+". Is the server working?" );
          response = {'exists':false,'isfile':false};
      }
    });
    console.log(response)
    return response;
}

/**
 *
 *
 *
 *
 **/
function create_folder(location){
    var response = undefined;
    $.ajax({
      url: "/create_directory",
      type: "POST",
      data: JSON.stringify({'location':location}),
      dataType: "text",
      async: false,
      complete: function(jqXHR, textStatus) {
          console.log(jqXHR.responseText)
          response = $.parseJSON(jqXHR.responseText);
      },
      
      error:function( jqXHR, textStatus, errorThrown ){
          alert( "Request failed: " + textStatus+". Is the server working?" );
          response = {'done':false};
      }
    });
    console.log(response)
    return response;
}
