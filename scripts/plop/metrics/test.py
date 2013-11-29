from metrics import processDir,mean,filterRecords,genSingleTraj,genMetricsFile

import sys

## extract files with energy lesser than -26770 && processor == 1

records = processDir(sys.argv[1],"traj")
selection = filterRecords("Proc == 1 and Energy<-26759",records)
genSingleTraj("file.pdb",records,selection)
genMetricsFile("metrics.dat",["Proc","Energy"],selection)