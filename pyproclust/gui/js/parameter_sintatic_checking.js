function check(this_rules, for_this_algorithms_parameters){
    for (var i =0; i< this_rules.length; i++){
        var params = [];
        var param_names = [];
        for (var param_key in this_rules[i]["params"]){
            param_names.push(param_key);
            params.push(this_rules[i]["params"][param_key]);
            console.log(param_key + " "+this_rules[i]["params"][param_key])
        }
        
        try{
            this_rules[i]["rule"](params);
        }
        catch(error_message){
            throw "Error in "+for_this_algorithms_parameters+" parameters ("+
                    param_names+"): "+error_message;
        }
    }
}

function lists_have_same_length(lists){
    var list1 = lists[0];
    var list2 = lists[1];
    
    if(list1.length != list2.length){
        throw "Lists must have the same length.";
    }
}

function list_not_empty(list){
    if(list[0].length == 0){
        throw "List cannot be empty."
    }
}
