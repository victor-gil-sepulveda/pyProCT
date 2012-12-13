'''
Created on 16/03/2012

@author: victor
'''
import subprocess
import prody #@UnresolvedImport

def getPDBStructure(pdb_path1,pdb_path2=None,selection_string="",subset_str = "calpha"):
    prody.setVerbosity('none')
    pdb = prody.parsePDB(pdb_path1, subset='calpha')
    if selection_string != "":
        selection = pdb.select(selection_string)
        pdb = selection.copy()
    if pdb_path2 != None:
        print "Found 2 pdbs, getting structures"
        pdb2 = getPDBStructure(pdb_path2,None,selection_string,subset_str)
        pdb.addCoordset(pdb2)
        print "Union has ",len(pdb.getCoordsets())," coordsets."
    return pdb

def get_number_of_frames(pdb_file):
    """
    Uses egrep to count the number of times MODEL appears in the trajectory file, which will be indeed
    the number of different structures.
    """
    process = subprocess.Popen(["egrep","-c",''"^MODEL"'',pdb_file],stdout=subprocess.PIPE)
    lines = process.stdout.readlines()
    return int(lines[0])

def read_to_TAG(file_in_handler,TAG):
    """
    Advances the file handler reading cursor to the line after and ENMDL or TER tags.
    It has unexpected behavior if both are found.
    """    
    lines = []  
    for l in file_in_handler:
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
    raise Exception, "No ENDMDL or TER tags found. This file may not be correct." 

def advance_to_TAG(file_in_handler,TAG):
    """
    Advances the file handler reading cursor until it finds a MODEL tag.
    """   
    for l in file_in_handler:
        if l[:len(TAG)]==TAG:
            return

def write_a_tfile_model_into_other_tfile(file_in_handler, file_out_handler, current_model, INITIAL_TAG, END_TAG, skip=False, keep_header = False):
    """
    It writes the next model found in a trajectory unless the skip parameter is set to 'True'.
    It makes the input file handler reading cursor to be placed the line after the end of 
    the model.
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
    

def extract_frames_from_trajectory(file_handler_in, number_of_frames, file_handler_out, frames_to_save, INITIAL_TAG="MODEL", END_TAG="ENDMDL", keep_header = False):
    """
    It extracts some frames from one trajectory and writes them down in an opened file handler.
    """
    current = 0
    for i in range(number_of_frames):
        if i in frames_to_save:
            write = True
            current = current + 1
        else:
            write = False
        write_a_tfile_model_into_other_tfile(file_handler_in,file_handler_out,current,INITIAL_TAG,END_TAG,not write, keep_header)

def create_CA_file(file_handler_in, file_handler_out):
    """
    Creates a valid CA pdb from a valid pdb. Useful to reduce load times.
    """
    for l in file_handler_in:
        if l[0:3] == 'MOD' or l[0:3] == 'END' or l[0:3] == 'TER': 
            file_handler_out.write(l)
        if l[0:3] == 'ATO':
            if l[12:16] ==" CA ":
                file_handler_out.write(l)
                
def repair_MODEL_ENDMDL_tags(input_pdb, output_pdb):
    """
    Tries to discover the structure of the file, to partition it in models. Then it adds the MODEL, ENDMDL tags
    for each of its discovered models, generating this way a correct PDB.
    
    @param input_pdb: Is the pdb file we want to repair.
    @param output_pdb: Is the resulting file (hopefully correct). 
    """
    file_in = open(input_pdb,"r")
    boundaries = []
    line_number = 0
    reading_model = False
    last_atom_serial = -1
    for line in file_in:
        if line[0:4] == "ATOM":
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
    file_in.close()
    
    file_in = open(input_pdb,"r")
    file_out = open(output_pdb,"w")
    line_number = 0
    current_boundary = 0
    for line in file_in:
        if current_boundary == len(boundaries):
            break
        #print current_boundary
        if line_number == boundaries[current_boundary][0]:
            file_out.write("MODEL "+str(current_boundary).rjust(8)+"\n")        
        if line_number>=boundaries[current_boundary][0] and line_number<=boundaries[current_boundary][1]:
            file_out.write(line)
        if line_number == boundaries[current_boundary][1]:
            file_out.write("ENDMDL\n")        
            current_boundary += 1
        line_number += 1
    
    file_in.close()
    file_out.close()