import sys

STRIDE = 6
END = 104
END_INCLUSIVE = END+1
for i in range(0,END_INCLUSIVE, STRIDE):
    file_range = (i,min(END_INCLUSIVE,i+STRIDE))
    print "Merging %d -> %d."%file_range
    output_handler = open("merged/%d_%d.pdb"%(file_range[0],file_range[1]-1),"w")
    for i in range(*file_range):
        file_handler = open("%d.pdb"%i,"r")
        for line in file_handler:
            output_handler.write(line)
        file_handler.close()
    output_handler.close()

