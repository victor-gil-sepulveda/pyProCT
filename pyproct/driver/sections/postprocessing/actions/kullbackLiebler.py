'''
Created on 16/07/2014

@author: victor
'''


class KullbackLieblerPostAction(object):
    KEYWORD = "kullback_liebler"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, trajectoryHandler, workspaceHandler, matrixHandler, generatedFiles):
        if KullbackLieblerPostAction.KEYWORD in postprocessing_parameters:
#             klDiv = KullbackLeiblerDivergence(self.trajectoryHandler.pdbs, self.matrixHandler.distance_matrix)
#             kl_file_path = os.path.join(self.workspaceHandler["matrix"], "kullback_liebler_divergence")
#             klDiv.save(kl_file_path)
#             matrix_image_file_path = os.path.join(self.workspaceHandler["matrix"], parameters["matrix"]["image"]["filename"])
#             self.generatedFiles.append({"description":"Kullback-Leibler divergence",
#                                         "path":matrix_image_file_path,
#                                         "type":"text"})