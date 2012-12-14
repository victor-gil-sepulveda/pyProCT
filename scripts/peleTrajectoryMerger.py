import os
import sys
import optparse
import subprocess

def get_number_of_frames(pdb_file):
    """
    Uses egrep to count the number of times MODEL appears in the trajectory file, which will be indeed
    the number of different structures.
    """
    process = subprocess.Popen(["egrep","-c",''"^MODEL"'',pdb_file],stdout=subprocess.PIPE)
    lines = process.stdout.readlines()
    return int(lines[0])

def write_a_model(file_handler,current_model,out_file_handler):
    """
    Precondition: All 'MODEL' tag has a 'ENDMDL' tag.
    """    
    # Find MODEL tag   
    l = file_handler.readline()
    while not l[:5]=="MODEL":
        l = file_handler.readline()
    
    # Once the model tag is found... write it down
    out_file_handler.write("MODEL"+str(current_model).rjust(9)+"\n")
        
    # Then write everything 'till the ENDMDL tag    
    l = file_handler.readline()
    while not l[:6]=="ENDMDL":
        out_file_handler.write(l)
        l = file_handler.readline()
    out_file_handler.write("ENDMDL\n")

def get_number_of_complete_models(pdb_file):
    process = subprocess.Popen(["egrep","-c",''"^ENDMDL"'',pdb_file],stdout=subprocess.PIPE)
    lines = process.stdout.readlines()
    return int(lines[0])
    
if __name__ == '__main__':
    
    parser = optparse.OptionParser(usage='%prog -p <arg> -o <arg> [-d <arg>]',
                                   version='1.0')
    parser.add_option('-p', action="store", dest = "traj_prefix",help="The prefix of the trajectory files to be merged.",metavar = "traj_")
    parser.add_option('-d', action="store", dest = "directory", default = "./", help="Directory where the trajectories are.",metavar = "./my_trajectories/")
    parser.add_option('-o', action="store", dest = "out_file",help="Name of the merged trajectory.",metavar = "trajectory.pdb")
    options, args = parser.parse_args()
    
    # Mandatory options
    if (not options.traj_prefix) or (not options.out_file):
        parser.error("Please specify the prefix (-p) and output trajectory file (-o).")
    
    # Pick files 
    traj_file = []
    files=os.listdir(options.directory) 
    for filename in files:
        if options.traj_prefix in filename:
            traj_file.append(options.directory+filename)
    file_out=open(options.out_file,'w')
    traj_file.sort()

    models_left = []
    total_files = len(traj_file)
    total_models = 0
    total_complete_models = 0
    pdb_file_handlers = []
    print "Discovering files"
    for pdb_filename in traj_file:
        number_of_models = get_number_of_frames(pdb_filename)
        number_of_complete_models = get_number_of_complete_models(pdb_filename)
        print "file:",pdb_filename,"frames:",number_of_models,"(%d)"%(number_of_complete_models)
        models_left.append(number_of_complete_models)
        total_models = total_models + number_of_models
        total_complete_models += number_of_complete_models
        pdb_file_handlers.append(open(pdb_filename,'r'))
    
    print "Total number of models:",total_models
    print "Total number of complete models:",total_complete_models

    current_model = 0
    # While there's any model left...
    while current_model < total_complete_models:
        # For each of the files write the next model
        for i in range(len(pdb_file_handlers)):
            if models_left[i] > 0:
                sys.stdout.flush()
                write_a_model(pdb_file_handlers[i],current_model, file_out)
                current_model = current_model + 1
                models_left[i] = models_left[i] - 1
                if current_model % 1000 == 0 :
                    print "Starting to process frame number", current_model
                    sys.stdout.flush()
    file_out.close()            