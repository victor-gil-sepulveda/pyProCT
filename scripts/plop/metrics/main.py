from metrics import processDir, genMetricsFile
import sys

records = processDir(sys.argv[1],"traj")
print records
genMetricsFile("metrics.dat", ["Energy","L1BindingEne"], records)