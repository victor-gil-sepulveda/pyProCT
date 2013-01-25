/*
    Shows and hides one field when the button is clicked, depending on its initial
    state.
*/
function toggable_checkbox(this_checkbox, over_this_field){
    $(this_checkbox).click(function() {
        $(over_this_field).toggle("blind")
    });
}

/* 
    This function adds a behaviour to a group of 2 radiobuttons.
    Each radiobutton has one related field which will be hidden or unhidden depending
    on it's radiobutton state. This means that when radiobutton 1 is checked, its field 
    shows up and the second radiobutton field is hidden (and viceversa).  
*/
function toggable_radiobutton(radio1, radio2, radio1_toggles_this, radio2_toggles_this){
    $(radio1).click(function() {
        do_toggling(radio1_toggles_this, radio2_toggles_this);
    });
    $(radio2).click(function() {
        do_toggling(radio2_toggles_this, radio1_toggles_this);
    });
}

/*
    Helper function.
    Given two fields which cannot be hidden at the same time, shows one and hides the
    other, updating the hidden parameter. Consecutive calls to this function with same 
    parameters doesn't do anything (as 'hidden' value is checked).
*/
function do_toggling(my_thing_to_show, the_thing_to_hide){
    if (my_thing_to_show != "" ){
        if($(my_thing_to_show).prop('hidden') == true){
            if (my_thing_to_show!= ""){
                $(my_thing_to_show).show({effect:"blind",
                    complete:function(){$(my_thing_to_show).prop('hidden',false)}});
            }
            if (the_thing_to_hide != ""){
                $(the_thing_to_hide).hide({effect:"blind",
                    complete:function(){$(the_thing_to_hide).prop('hidden',true)}});
            }
        }
    }
    else{
        if (the_thing_to_hide != ""){
            $(the_thing_to_hide).hide({effect:"blind",
                complete:function(){$(the_thing_to_hide).prop('hidden',true)}});
        }
    }
}
