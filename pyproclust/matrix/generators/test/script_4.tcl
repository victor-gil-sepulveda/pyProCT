set inpfile "a_trajectory.pdb"

set outfile  "partial_out_4"

set outdescriptor  [open "${outfile}" w]

set starting_frame 4

set total_frames  1

set sym [lindex $argv 4]

set molid [mol load pdb $inpfile] 

#set nf [molinfo $molid get numframes]

set nf 8

puts "Calculating the RMSD matrix for $nf frames..."
set last_frame [expr $total_frames + $starting_frame]
for { set f1 $starting_frame } { $f1 < $last_frame } { incr f1 1 } {
    
    puts "Calculating frame $f1 of $nf"        
    
    set selection_f1 [atomselect $molid "name CA" frame $f1]
    set fit_selection_f1 [atomselect $molid "name CA" frame $f1]
    set time_start [clock clicks -milliseconds]
    # Set second start to i+1    
    set second_start $f1
    incr second_start 1  
    for { set f2 $second_start } { $f2 < $nf } { incr f2 1 } {
        set selection_f2 [atomselect $molid "name CA" frame $f2]    
        set fit_selection_f2 [atomselect $molid "name CA" frame $f2]
        
        $fit_selection_f2 move [measure fit $selection_f2 $selection_f1]

        puts -nonewline $outdescriptor [format "%.4f " [measure rmsd $fit_selection_f2 $selection_f1]]

        $selection_f2 delete
        $fit_selection_f2 delete
    }
    puts $outdescriptor ""
    
    puts "It took [expr [clock clicks -milliseconds] - $time_start] ms."

    $selection_f1 delete
    $fit_selection_f1 delete
}
puts "Script finished. $total_frames processed. Closing vmd..."
close $outdescriptor
exit
