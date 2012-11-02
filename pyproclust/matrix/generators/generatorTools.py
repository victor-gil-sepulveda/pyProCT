import math

def get_symm_ranges(number_of_processors, number_of_frames):
    """
    Generates some ranges so that for a given number of processors, each one processes
    as many operations as possible while keeping the balance.
    A range contains two values. In the rmsd matrix context the first value represents
    the number of frame of the trajectory. The second would be the amount of frames
    to be processed. This ranges are calculated to balance the work of creating a square 
    symmetric matrix.
    """
    ranges = []
    
    starting_frame = 0
    frames_per_iteration = int(float(number_of_frames) / number_of_processors)
    
    # some rounding
    if float(number_of_frames) / number_of_processors - frames_per_iteration > 0.5 :
        frames_per_iteration = frames_per_iteration + 1
    
    
    frames_left = number_of_frames
    
    end_condition = True
    
    while(end_condition): 
        
        frames_this_iteration =  frames_per_iteration
        
        if frames_left  > frames_per_iteration and \
        frames_left  < 2*frames_per_iteration :
            frames_this_iteration =  number_of_frames - starting_frame
        
        if frames_this_iteration > 0:
            ranges.append((starting_frame,frames_this_iteration))
        
        starting_frame = starting_frame + frames_per_iteration
        frames_left = frames_left - frames_this_iteration
    
        end_condition = frames_left > 0
    
    return ranges


def get_nosymm_ranges(number_of_processors, number_of_frames):
    """
    It does exactly the same that  get_symm_ranges, but in this case the work is shared
    to make all processors process a similar number of elements. That's because in this 
    case the goal is generate the upper triangle of a symmetric matrix.
    """
    ranges = []
    
    num_elements = number_of_frames*(number_of_frames-1)/2.
    
    num_elements_per_processor = int(math.floor(float(num_elements) / number_of_processors))
    
    elements_this_line = number_of_frames - 1 
    
    starting_frame  = 0 
    for i in range(number_of_processors): #@UnusedVariable
        elements_for_this_processor = 0
        lines = 0
        
        while (elements_for_this_processor < num_elements_per_processor and elements_this_line > 0): 
            elements_for_this_processor += elements_this_line
            elements_this_line = elements_this_line - 1 
            lines += 1
            
        ranges.append((starting_frame,lines))
        starting_frame = starting_frame + lines
        
    return ranges