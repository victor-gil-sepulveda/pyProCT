import sys

file_range = (int(sys.argv[1]),int(sys.argv[2])+1)
output_handler = open(sys.argv[3],"w")
for i in range(*range):
    file_handler = open(str(i),"r")
    for line in file_handler:
        output_handler.write(line)
    file_handler.close()
output_handler.close()

