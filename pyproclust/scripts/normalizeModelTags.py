import sys

file_in = open(sys.argv[1],"r")
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


file_out = open(sys.argv[2],"w")
file_in = open(sys.argv[1],"r")
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
