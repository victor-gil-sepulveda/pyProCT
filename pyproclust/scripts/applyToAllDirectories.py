import os
import random
from pyproclust.tools.pdbTools import get_number_of_frames, extract_frames_from_trajectory


directories = []
for dirname, dirnames, filenames in os.walk('.',topdown=False):
    try:
        directories.append(dirname.split("/")[1])
    except:
        pass    
directories =  list(set(directories))

types = ['_5_0.75','_5_1.4','_10_0.75','_10_1.4']
#
#script = '/home/victor/workspaces/Python/TrajectoryClustering/src/pyproclust/scripts/peleTrajectoryMerger.py -p traj_ -o %s_merged.pdb -d'
#for d in directories:
#    print "Processing ", d
#    os.system("pwd")    
#    print "python "+script%(d)+" "+d+'/'
#    os.system("python "+script%(d)+" "+d+'/')

# create 10k frames mixed trajectories (33% from pele_amber [pele with a starting frame from amber]
# , 33" from pele_charmm and 33% from pele_opls)

#for t in types:
#    print "--------------------"
#    print "Type ", t
#    print "--------------------"
#    selected = []
#    for f in filenames:
#        if t in f:
#            selected.append(f)
#    filename = "merged_extracted"+t+".pdb"
#    file_out = open(filename,'w')
#    for s in selected:
#        print "Working in ", s
#        num_frames_to_extract = 3334
#        num_frames_in_file = get_number_of_frames(s)
#        frame_numbers = range(num_frames_in_file)
#        random.shuffle(frame_numbers)
#        frame_extraction = frame_numbers[:num_frames_to_extract]
#        file_in = open(s,'r')
#        extract_frames_from_trajectory(file_in,num_frames_in_file,file_out,frame_extraction)
#        file_in.close()
#    file_out.close()
#    print "A pdb file with ",get_number_of_frames(filename), " frames was created."

## Create files of at most 10k elements of each (10ks of pele with amber as first frame, etc...)
#first_frame_types  = ["amber","charmm","opls"]
#for t in first_frame_types:
#    print "--------------------"
#    print "Type ", t
#    print "--------------------"
#    selected = []
#    for f in filenames:
#        if t in f:
#            selected.append(f)
#    filename = t+"_10ks.pdb"
#    file_out = open(filename,'w')
#    for s in selected:
#        print "Working in ", s
#        num_frames_to_extract = 10000 / len(selected)
#        num_frames_in_file = get_number_of_frames(s)
#        frame_numbers = range(num_frames_in_file)
#        random.shuffle(frame_numbers)
#        frame_extraction = frame_numbers[:num_frames_to_extract]
#        file_in = open(s,'r')
#        extract_frames_from_trajectory(file_in,num_frames_in_file,file_out,frame_extraction)
#        file_in.close()
#    file_out.close()
#    print "A pdb file with ",get_number_of_frames(filename), " frames was created."
#
#exit()

base_script = """
from pyproclust.protocol.protocolImplementation import Protocol
from pyproclust.protocol.protocolParameters import ProtocolParameters


if __name__ == '__main__':
    
    # Example
    proto_params =  ProtocolParameters()
    
    # Some common stuff
    proto_params.number_of_processors = 4
    proto_params.pdb1 = "%s"
    proto_params.pdb2 = "%s"
    proto_params.rmsd_selection = "resnum < 71"
    
    # Matrix Creation / Loading
    proto_params.store_matrix_path = "matrix"
    #proto_params.matrix_file = "matrix"
    
    # Clustering algorithms parameters
    proto_params.use_kmedoids = True
    proto_params.use_hierarchical = True
    proto_params.use_dbscan = True
    proto_params.hierarchical_cutoff_list = [] # Please calculate this cutoffs for me
    proto_params.kmedoids_step = 10
    proto_params.dbscan_param_pairs = [] # Please search them for me
    
    # Refinement
    proto_params.do_refinement = True
    proto_params.refinement_min_clusters = 10
    proto_params.refinement_max_clusters = 60
    proto_params.refinement_step = 5
    
    # Cluster selection
    proto_params.max_noise = 15
    proto_params.min_clusters = 10
    proto_params.max_clusters = 30
    proto_params.min_cluster_size = 50
    proto_params.evaluation_types = [ "NumClusters","NumClusteredElems","MeanClusterSize","NoiseLevel","CythonMirrorCohesion",\
                "CythonSilhouette","PCAanalysis","NormNCut"]
    proto_params.cluster_score_value_map = [{"CythonSilhouette":(0.9,">"),"PCAanalysis":(1.,"<"),"NormNCut":(0.3,">")},\
                               {"CythonSilhouette":(0.8,">"),"PCAanalysis":(1.,"<"),"CythonMirrorCohesion":(0.2,"<"),"NormNCut":(0.1,">")}]
    # Results params
    proto_params.report_file = "report"
    proto_params.comparator_results_path = "overlap.txt"
   
    ####################
    # Start partying
    ####################
    protocol = Protocol()
    protocol.run(proto_params)
"""

trajectories = {'amber':"/home/victor/Escritorio/TestDatasets/MDs/amber_corrected.pdb",
                'charmm':"/home/victor/Escritorio/TestDatasets/MDs/charmm.pdb",
                'opls':"/home/victor/Escritorio/TestDatasets/MDs/opls.pdb",
                'pele_amber':"/home/victor/Escritorio/TestDatasets/PELE/amber_10ks.pdb",
                'pele_charmm':"/home/victor/Escritorio/TestDatasets/PELE/charmm_10ks.pdb",
                'pele_opls':"/home/victor/Escritorio/TestDatasets/PELE/opls_10ks.pdb",
                'mixed_pele_5_0.75':"/home/victor/Escritorio/TestDatasets/PELE/merged_extracted_5_0.75.pdb",
                'mixed_pele_5_1.4':"/home/victor/Escritorio/TestDatasets/PELE/merged_extracted_5_1.4.pdb",
                'mixed_pele_10_0.75':"/home/victor/Escritorio/TestDatasets/PELE/merged_extracted_10_0.75.pdb",
                'mixed_pele_10_1.4':"/home/victor/Escritorio/TestDatasets/PELE/merged_extracted_10_1.4.pdb"}

trajectory_types = trajectories.keys()
                # Comprobar si el frame inicial sesga la simulacion
comparisons = [ ('amber','pele_amber'),('amber','pele_charmm'),('amber','pele_opls'),\
                ('charmm','pele_amber'),('charmm','pele_charmm'),('charmm','pele_opls'),\
                ('opls','pele_amber'),('opls','pele_charmm'),('opls','pele_opls'),                                                                                
                # Comprobar la influencia de los parametros                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                ('amber','mixed_pele_5_0.75'),('amber','mixed_pele_5_1.4'),('amber','mixed_pele_10_0.75'),('amber','mixed_pele_10_1.4'),
                ('charmm','mixed_pele_5_0.75'),('charmm','mixed_pele_5_1.4'),('charmm','mixed_pele_10_0.75'),('charmm','mixed_pele_10_1.4'),
                ('opls','mixed_pele_5_0.75'),('opls','mixed_pele_5_1.4'),('opls','mixed_pele_10_0.75'),('opls','mixed_pele_10_1.4')]

#MD
comparisons = [('amber','charmm'), ('amber','opls'), ('charmm','opls'), ("opls","mixed_pele_5_0.75"), ("charmm","mixed_pele_10_1.4")]

for type1,type2 in comparisons:
    folder = type1+"_vs_"+type2
    os.system("mkdir "+folder)
    script_file = open(folder+"/script.py","w")
    script_file.write(base_script%(trajectories[type1],trajectories[type2]))                                                                
    script_file.close()
    os.system("rm -rf "+folder+"/clusterings")
    os.system("rm -rf "+folder+"/refinement")
    os.chdir(folder)
    os.system("python script.py")                                                         
    os.chdir("..")
        
#for traj in trajectories:
#    for typ in types:
#        folder = "comp_"+traj+"vspele"+typ
#        os.system("mkdir "+folder)
#        pele_filename = "merged_and_extracted"+typ+".pdb"
#        script_file = open(folder+"/script.py","w")
#        script_file.write(base_script%(trajectories[traj],pele_filename))
#        script_file.close()
#        os.system("cp MD/"+trajectories[traj]+" "+folder)
#        os.system("cp "+pele_filename+" "+folder)
