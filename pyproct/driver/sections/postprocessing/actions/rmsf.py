'''
Created on 16/07/2014

@author: victor
'''
from pyproct.tools import visualizationTools
import json
import os


class RmsfPostAction(object):
    KEYWORD = "rmsf"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):
        if RmsfPostAction.KEYWORD in postprocessing_parameters:
#             try:
                displacements_path, CA_mean_square_displacements = visualizationTools.calculate_RMSF(clustering,
                                                                                                     trajectoryHandler,
                                                                                                     workspaceHandler,
                                                                                                     matrixHandler)

                generatedFiles.append({
                                            "description":"Alpha Carbon mean square displacements",
                                            "path":os.path.abspath(displacements_path),
                                            "type":"text"
                })

                open(displacements_path,"w").write(json.dumps(CA_mean_square_displacements,
                                                      sort_keys=False,
                                                      indent=4,
                                                      separators=(',', ': ')))
#             except Exception:
#                 print "[ERROR][Driver::postprocess] Impossible to calculate rmsf file."