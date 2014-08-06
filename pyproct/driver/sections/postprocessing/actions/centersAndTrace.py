"""
Created on 16/07/2014

@author: victor
"""
from pyproct.tools import visualizationTools
import json
import os

class CentersAndTracePostAction(object):
    KEYWORD = "centers_and_trace"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):
        if CentersAndTracePostAction.KEYWORD in postprocessing_parameters:
            try:
                centers_path, centers_contents = visualizationTools.generate_selection_centers_file(clustering,
                                                                                                    workspaceHandler,
                                                                                                    trajectoryHandler)

                open(centers_path,"w").write(json.dumps(centers_contents,
                                          sort_keys=False,
                                          indent=4,
                                          separators=(',', ': ')))
                generatedFiles.append({
                                            "description":"Centers of the selection used to calculate distances",
                                            "path":os.path.abspath(centers_path),
                                            "type":"text"
                })

            except Exception:
                print "[ERROR][Driver::postprocess] Impossible to calculate selection centers file."