import os.path
from pyproct.tools.commonTools import convert_to_utf8
import json

redo_compression_script = """
{
     "global": {
        "control": {
            "scheduler_type": "Serial"
        },
        "workspace": {
            "base": "redo_clust_tmp"
        }
    },
    "data": {
        "files": [
            {
                "atoms_file": "/gpfs/scratch/bsc72/bsc72476/Victor/2JOF/2J0F.pdb",
                "file": "%s"
            }
        ],
        "type": "pdb_ensemble",
        "matrix": {
            "method": "rmsd",
            "parameters": {
                "calculator_type": "QCP_OMP_CALCULATOR",
                "fit_selection": "name CA"
            }
        }
    },
    "clustering": {
        "generation": {
            "method": "load",
            "clusters":%s
        }
    },
    "postprocess": {
        "compression": {
            "final_number_of_frames": 2000,
            "type": "KMEDOIDS"
        }
    }
}

"""
def load_dic_in_json(filename):
    return convert_to_utf8(json.loads("".join(open(filename,"r").readlines())))

def get_best_clustering(results_file):
    results = load_dic_in_json(results_file)
    best_clustering_id = results['best_clustering']
    return results['selected'][best_clustering_id]

all_initial_files = open("files.txt")
BASE_TRAJ_FOLDER = "/gpfs/scratch/bsc72/bsc72476/Victor/2JOF"

i = 0
for file_path in all_initial_files.readlines():
    traj_path = os.path.join(BASE_TRAJ_FOLDER, file_path.strip())
    clustering_results_path = os.path.join("compressions","level2","clustering_%d"%i,"results","results.json")
    clusters = get_best_clustering(clustering_results_path)["clustering"]["clusters"]
    script =  (redo_compression_script%(traj_path,str(clusters))).replace("'",'"')
    open("tmp_script.json","w").write(script)
    os.system("python -m pyproct.main tmp_script.json")
    os.system("cp redo_clust_tmp/results/compressed_pdb.pdb redone_comp/%d.pdb"%i)
    i = i+1
all_initial_files.close()