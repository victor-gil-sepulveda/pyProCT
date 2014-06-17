import pyproct.tools.scriptTools as tools
import os.path
import sys

base_path = sys.argv[1]
file_paths = sorted(os.listdir(base_path))
# ["compressions/level2/clustering_0",
# "compressions/level2/clustering_1",
# "compressions/level2/clustering_10",
# "compressions/level2/clustering_11",
# "compressions/level2/clustering_12",
# "compressions/level2/clustering_13",
# "compressions/level2/clustering_14",
# "compressions/level2/clustering_15",
# "compressions/level2/clustering_16",
# "compressions/level2/clustering_17",
# "compressions/level2/clustering_18",
# "compressions/level2/clustering_19",
# "compressions/level2/clustering_2",
# "compressions/level2/clustering_20",
# "compressions/level2/clustering_21",
# "compressions/level2/clustering_22",
# "compressions/level2/clustering_23",
# "compressions/level2/clustering_24",
# "compressions/level2/clustering_26",
# "compressions/level2/clustering_27",
# "compressions/level2/clustering_28",
# "compressions/level2/clustering_29",
# "compressions/level2/clustering_3",
# "compressions/level2/clustering_30",
# "compressions/level2/clustering_31",
# "compressions/level2/clustering_32",
# "compressions/level2/clustering_34",
# "compressions/level2/clustering_35",
# "compressions/level2/clustering_36",
# "compressions/level2/clustering_37",
# "compressions/level2/clustering_4",
# "compressions/level2/clustering_42",
# "compressions/level2/clustering_43",
# "compressions/level2/clustering_44",
# "compressions/level2/clustering_45",
# "compressions/level2/clustering_46",
# "compressions/level2/clustering_47",
# "compressions/level2/clustering_48",
# "compressions/level2/clustering_5",
# "compressions/level2/clustering_6",
# "compressions/level2/clustering_7",
# "compressions/level2/clustering_8",
# "compressions/level2/clustering_9"]

total_done = 0
for path in file_paths:
    complete = os.path.join(base_path, path,"results", "compressed_pdb.pdb")
    if not os.path.exists(complete):
        print complete
    else:
        total_done += 1

print "completed ",total_done