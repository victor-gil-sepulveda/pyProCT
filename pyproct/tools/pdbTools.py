"""
Created on 16/03/2012

@author: victor
"""

def get_model_tags(pdb_file):
    """
    Uses the 'egrep' shell command to count the number of times MODEL appears in the trajectory file, which will be indeed
    the number of different structures.
    
    :param pdb_file: The pdb file (path) from which the number of frames will be counted.
    
    :return: The lines containing the found model tags.
    """
    model_lines = []
    for line in open(pdb_file, "r"):
        if line[0:5] == "MODEL":
            model_lines.append(line)
    return model_lines

def get_number_of_frames(pdb_file):
    """
    Uses the 'egrep' shell command to count the number of times MODEL appears in the trajectory file, which will be indeed
    the number of different structures.
    
    :param pdb_file: The pdb file (path) from which the number of frames will be counted.
    
    :return: The number of MODEL tags found.
    """
    return len(get_model_tags(pdb_file))

def get_remarks(pdb_file):
    """
    Sometimes pdb files have remarks with extra information before the MODEL tag. 
    This remarks have a fixed format, however sometimes, when adding user information, this 
    format is not honored.
    :param pdb_file: The pdb file (path) from which the number of frames will be counted.
    
    :return: an array of arrays with the remarks (one array of lines per model). 
    
    """
  
    handler = open(pdb_file, "r")
    line_groups = []
    lines = []
    starting_atom = True
    
    for line in handler:
        if line[0:6] == "REMARK":
            lines.append(line)
            starting_atom = True
        
        if line[0:5] == "MODEL" or (line[0:4] == "ATOM" and starting_atom):
            line_groups.append(lines)
            lines = []
            starting_atom = False
    
    handler.close()
    
    return line_groups

def filter_remarks(remarks_list, subset="NONE"):
    """
    Filters a remarks list using a subset criteria.
    :param remarks_list: the list we want to filter.
    :param subset: The criteria to do the filtering. Possible values are:
            - "NONE": not to store remarks (Default)
            - "STANDARD": stores remarks that follow pdb standard
            - "NOT STANDARD": stores remarks not following the pdb standard
            - "ALL": stores all remarks
    :return: The filtered list of remarks. 
    """
    available =  ["NONE","ALL","STANDARD", "NOT STANDARD"]
    if not subset in available:
        print "[WARNING proteinEnsembleDataLoader::filter_remarks] %s subset is not an available option %s."%(subset, str(available))
        return []
    
    filtered_groups = []
    for group in remarks_list:
        filtered_remarks = []
        for remark in group:
            standard_remark = (remark[6] == " " and remark[7:10].isdigit() and remark[10] == " ")
            if subset == "ALL":
                filtered_remarks.append(remark)
            
            elif subset == "STANDARD" and standard_remark: 
                filtered_remarks.append(remark)
            
            elif subset == "NOT STANDARD" and not standard_remark:
                filtered_remarks.append(remark)
        filtered_groups.append(filtered_remarks)
    return filtered_groups

def get_number_of_atoms(pdb_file):
    """
    Counts the number of ATOM and HETATM entries in the first model of a pdb file.
    
    @param pdb_file: The pdb file (path) from which the number of atoms will be counted.
    
    @return: The number of atoms found in the first frame.
    """
    handler = open(pdb_file,"r")
    number_of_models = 0
    number_of_atoms = 0
    for line in handler:
        
        if line[0:5] == "MODEL":
            number_of_models += 1
            if number_of_models > 1:
                break

        if number_of_models == 1 and (line[0:4] == "ATOM" or line[0:6] == "HETATM"):
            number_of_atoms += 1
    return number_of_atoms

def read_to_TAG(file_input_handler, TAG):
    """
    Advances the file handler reading cursor to the line after a tag.
    
    @param file_input_handler: Working file handler.
    @param TAG: The tag to advance to (for instance "MODEL", or "TER"). 
    """
    lines = []  
    for l in file_input_handler:
        if isinstance(TAG, basestring):
            if l[:len(TAG)]==TAG:
                return lines
            else:
                lines.append(l)
        else: # if not, is a list
            for tag in TAG:
                if l[:len(tag)]==TAG:
                    return lines
                lines.append(l)
    raise Exception, "The end of the file was reached and the tag "+TAG+" was not found." 

def advance_to_TAG(file_input_handler, TAG):
    """
    Advances the file handler reading cursor until it finds a tag.
    
    @param file_input_handler: Working file handler.
    @param TAG: The tag to advance to (for instance "MODEL", or "TER"). 
    """   
    for l in file_input_handler:
        if l[:len(TAG)]==TAG:
            return

def write_a_tfile_model_into_other_tfile(file_in_handler, 
                                         file_out_handler, 
                                         current_model, 
                                         INITIAL_TAG, 
                                         END_TAG, 
                                         skip=False, 
                                         keep_header = False):
    """
    It writes the next model found in a trajectory unless the skip parameter is set to 'True'.
    It makes the input file handler reading cursor to be placed the line after the end of the model.
    
    @param file_in_handler: Input file handler.
    @param file_out_handler: File handler to write the model to.
    @param current_model: Is the model number we will write in the output handler.
    @param INITIAL_TAG: Is the tag preceding the atom pdb data. Usually is "MODEL".
    @param END_TAG: Is the last tag in a model, after the pdb atom data. Usually is "ENDMDL" or "TER".
    @param skip: If true, the model will not be written to the output file, but the cursor of the input_file will be moved
    anyway.
    @param keep_header: Will try to keep any header previous to INITIAL_TAG, such as "REMARK" lines.
    """    
    # Shall we write this model or not?
    if skip:
        read_to_TAG(file_in_handler,END_TAG)
    else: 
        if not keep_header:
            # Find model tag
            advance_to_TAG(file_in_handler,INITIAL_TAG)
        else:
            header = read_to_TAG(file_in_handler,INITIAL_TAG)
            file_out_handler.writelines(header)
        
        file_out_handler.write("MODEL"+str(current_model).rjust(9)+"\n")
            
        # Once the model tag is found... write it down
        lines  = read_to_TAG(file_in_handler,END_TAG)
        file_out_handler.writelines(lines)
        file_out_handler.write("ENDMDL\n")
    

def extract_frames_from_trajectory_sequentially(file_handler_in, 
                                   number_of_frames, 
                                   file_handler_out, 
                                   frames_to_save, 
                                   INITIAL_TAG="MODEL", 
                                   END_TAG="ENDMDL", 
                                   keep_header = False,
                                   write_frame_number_instead_of_correlative_model_number = False):
    """
    It extracts some frames from one trajectory and writes them down in an opened file handler.
    
    @param file_in_handler: Input file handler.
    @param number_of_frames: Number of models in the input trajectory.
    @param file_out_handler: File handler to write the model to.
    @param frames_to_save: An array with the frame numbers to copy in the output trajectory. For example: [0,1,2] will store the
    three first frames of the input trajectory into the second.
    @param INITIAL_TAG: Is the tag preceding the atom pdb data. Usually is "MODEL".
    @param END_TAG: Is the last tag in a model, after the pdb atom data. Usually is "ENDMDL" or "TER".
    @param keep_header: Will try to keep any header previous to INITIAL_TAG, such as "REMARK" lines.
    @param use_frame_number_as_model: If true, it will use the number of frame in 'frames_to_save' as model number
    in the new file. If not, models of the new file will have sequential numbering.
    """
    current = 0
    for i in range(number_of_frames):
        if i in frames_to_save:
            write = True
            if write_frame_number_instead_of_correlative_model_number:
                current = i
            else:
                current = current + 1
        else:
            write = False
        write_a_tfile_model_into_other_tfile(file_handler_in,
                                             file_handler_out,
                                             current,
                                             INITIAL_TAG,
                                             END_TAG,
                                             not write,
                                             keep_header)

def create_CA_file(file_handler_in, file_handler_out):
    """
    Creates a valid CA pdb from a valid pdb. Useful to reduce load times.
    
    @param file_handler_in: A file handler containing the trajectory to 'compress'.
    @param file_handler_out: The file handler in which the new trajectory is written.
    """
    for l in file_handler_in:
        if l[0:3] == 'MOD' or l[0:3] == 'END' or l[0:3] == 'TER': 
            file_handler_out.write(l)
        if l[0:3] == 'ATO':
            if l[12:16] ==" CA ":
                file_handler_out.write(l)
                
def get_model_boundaries(input_pdb_handler):
    """
    Tries to discover the structure of the file.
    
    @param input_pdb: Is the pdb file we want to repair.
    
    @return: The list of sections of the file that represent a model.
    """
    boundaries = []
    line_number = 0
    reading_model = False
    last_atom_serial = -1
    for line in input_pdb_handler:
        if line[0:4] == "ATOM" or line[0:4] == "HETA":
            atom_serial = int(line[5:11])
            # if the atom serial number decreases suddenly, it means we are starting
            # to read another frame of the trajectory, so we have to finish the last
            # and open a new boundary
            # if reading model was false, we are starting a new frame of the trajectory
            if last_atom_serial > atom_serial or reading_model == False:
                boundaries.append([line_number])
                reading_model = True
           
            last_atom_serial = atom_serial
        else:
            # if the tag is not atom, we are not reading a frame,
            # and we can close a boundary
            if reading_model == True:
                boundaries[-1].append(line_number-1)
                reading_model = False 
        
        line_number += 1
    return boundaries

def repair_MODEL_ENDMDL_tags(input_pdb_handler, output_pdb_handler, boundaries):
    """
    Uses the previously discovered partition with 'get_model_boundaries' to add the necessary MODEL, 
    ENDMDL tags for each of its discovered models.
    
    @param input_pdb: Is the pdb file handler we want to repair.
    @param output_pdb: Is the resulting file handler (hopefully correct). 
    """  
    line_number = 0
    current_boundary = 0
    for line in input_pdb_handler:
        if current_boundary == len(boundaries):
            break
        #print current_boundary
        if line_number == boundaries[current_boundary][0]:
            output_pdb_handler.write("MODEL "+str(current_boundary).rjust(8)+"\n")        
        if line_number>=boundaries[current_boundary][0] and line_number<=boundaries[current_boundary][1]:
            output_pdb_handler.write(line)
        if line_number == boundaries[current_boundary][1]:
            output_pdb_handler.write("ENDMDL\n")        
            current_boundary += 1
        line_number += 1

def grab_existing_frame_from_trajectory(trajectory_file_handler, output_file_handler, model_number):
    """
    Extracts a model of a trajectory and writes it into another file. The model must exist in the trajectory,
    otherwise the behavior will be undefined.
    
    @param trajectory_file_handler: Opened file handler containing the whole trajectory.
    
    @param output_file_handler: Opened file handler to write resulting model.
    
    @param model_number: The model we want to extract from the trajectory. It MUST exist.
    """
    write_model = False
    
    for line in trajectory_file_handler:
        if line[0:5] == "MODEL":
            if int(line[5:]) == model_number:
                write_model = True
        
        if write_model == True:
            output_file_handler.write(line)
        
        if write_model == True and  line[0:6] == "ENDMDL":
            write_model = False
    