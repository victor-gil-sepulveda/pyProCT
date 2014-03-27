'''
Created on 25/03/2014

@author: victor
'''
import os

PELE_MERGER = '/home/victor/workspaces/Python/pyProClust/scripts/pele/peleTrajectoryMerger.py'

def get_files(directory):
    files=os.listdir(directory)
    return sorted(files)

def gen_pdb_folders(to_this_temp):
    pdb_folders = []
    for i in range(to_this_temp+1):
        number = "{:0>2d}".format(i)
        pdb_folders.append("t%s/pdbs"%number)
    return pdb_folders

def gen_fixed_pdb_folders(to_this_temp):
    pdb_folders = []
    for i in range(to_this_temp+1):
        number = "{:0>2d}".format(i)
        pdb_folders.append("t%s_fixed"%number)
    return pdb_folders

if __name__ == '__main__':
    for folder, fixed_folder in zip(gen_pdb_folders(15),gen_fixed_pdb_folders(15)):
        os.system("mkdir %s"%fixed_folder)

        print "Working with %s ..."%folder

        for file_num, file_name in enumerate(get_files(folder)):
            fixed_file = "%s/%s.pdb"%(fixed_folder, "{:0>6d}".format(file_num))
            os.system("printf \"MODEL        0\n\" > %s" %(fixed_file))
            os.system("egrep \"^ATOM\" %s >> %s"%(folder+"/"+file_name, fixed_file))
            os.system("printf \"ENDMDL\n\" >> %s" %(fixed_file))
            file_num += 1
        print "\t-Fixed files created."

        MAX_FILES = 1000
        parsed_files = 0
        index_file = open("index.txt","w")
        index_files = ["index.txt"]
        while parsed_files < len(get_files(folder)):
            index_file.write(fixed_folder+"/{:0>6d}.pdb\n".format(parsed_files))
            parsed_files += 1
            if parsed_files%MAX_FILES == 0:
                index_file.close()
                file_name ="index_%s_%d.txt"%(fixed_folder, parsed_files/MAX_FILES)
                index_file = open(file_name,"w")
                index_files.append(file_name)
        index_file.close()

        print "\t-Index file created"
        for index_file in index_files:
            os.system("python %s -o %s.pdb -f %s --merge_action SEQUENTIAL --append"%(PELE_MERGER, fixed_folder, index_file))

        print "\t-File merged"

