function get_eval_dialog_contents(criteria_list){
    
    var contents = "\
        <table>\
        <tr><th>Action</th><th>Criteria</th><th>Weight</th></tr>\
        <tbody id = 'criteria_table'>"+
         get_table_contents_for_all_criteria(criteria_list)+
        "</tbody>\
        </table>";
    
    return contents;
 }
 
 function get_table_contents_for_all_criteria(criteria_list){
    var total_len = criteria_list.length;
    var contents = "";
    for (var i =0; i<total_len; i++){
        var row = get_row_contents(criteria_list[i]);
        contents += wrap_with_tr(row);
    }
    return contents;
 }
 
 function wrap_with_tr(this_html){
    contents = "<tr>"+this_html+"</tr>";
    return contents;
 }
 
 function get_row_contents(criteria_name){
    contents =  //wrap_with_td("<input class='dialog_checkbox' type='checkbox'/>", criteria_name)+
                wrap_with_td("<select>\
                                 <option> Minimize </option>\
                                 <option> Maximize </option>\
                            </select>", criteria_name)+
                wrap_with_td(criteria_name) +
                wrap_with_td("<input class='dialog_spinner "+criteria_name+"' type = 'input' value ='0'/>", criteria_name);
    return contents;
 }
 
 function wrap_with_td(this_html, this_class){
    contents = "<td class='"+this_class+"'>"+this_html+"</td>";
    
    return contents;
 }
 
function criteria_to_string(dialog_to_extract_data, criteria_list){
    var string_criteria = "";
    for (var i =0; i<criteria_list.length; i++){
        var weight =  $("."+criteria_list[i]+" .dialog_spinner").spinner( "value" );
        if (weight != 0){
            // Minimize or minimize?
            var min_max = $("."+criteria_list[i]+" option:selected").text();
            string_criteria += min_max+" "+criteria_list[i]+" (weigth: "+weight+") and"   
        }
    }
    // remove last and return
    return string_criteria.substring(0,string_criteria.length-4);
}

/*function criteria_to_string(dialog_to_extract_data, criteria_list){
    var string_criteria = "";
    for (var i =0; i<criteria_list.length; i++){
        var checked = $("."+criteria_list[i]+" input:checked");
        if (checked.length != 0){
            // Minimize or minimize?
            var min_max = $("."+criteria_list[i]+" option:selected").text();
            var weight =  $("."+criteria_list[i]+" .dialog_spinner").spinner( "value" );
            string_criteria += min_max+" "+criteria_list[i]+" (weigth: "+weight+") and"   
        }
    }
    // remove last and return
    return string_criteria.substring(0,string_criteria.length-4);
}*/
