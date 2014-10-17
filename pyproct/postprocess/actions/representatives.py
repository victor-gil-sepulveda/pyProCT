"""
Created on 10/02/2014

@author: victor
"""
import os.path
import prody
import cStringIO
from pyproct.tools.pdbTools import filter_remarks

class RepresentativesPostAction(object):

    KEYWORD = "representatives"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, data_handler, workspaceHandler, matrixHandler, generatedFiles):

        medoids = clustering.get_medoids(matrixHandler.distance_matrix)
        
        ids = [cluster.id for cluster in clustering.clusters]
    
        pdb_name = postprocessing_parameters.get_value("filename", default_value = "representatives")
        
        representatives_file_path = os.path.join( workspaceHandler["results"],"%s.pdb"%pdb_name)
    
        save_cluster_elements(medoids,
                              ids, 
                              representatives_file_path,
                              data_handler,
                              postprocessing_parameters)
    
        generatedFiles.append({
                               "description":"Cluster representatives",
                               "path":os.path.abspath(representatives_file_path),
                               "type":"pdb"
        })

def save_cluster_elements(elements,
                          ids,
                          out_pdb_name,
                          data_handler,
                          options):
    """
    Saves a pdb file containing the most representative elements of the clustering.

    @param elements: A list of the representative elements of the clustering we want to extract.

    @params ids: A list with the cluster ids (1 to 1 mapping with 'elements').

    @param out_pdb_name: The complete path of the produced file.

    @param data_handler: The trajectory handler for this run or an array with pdb file paths.

    @param options: postprocessing options to generate the file. Currently a dic with any of these:
        "keep_remarks" - Will add each model's remarks before the model header if present
            Possible values are:
            - "NONE": not to store remarks (Default)
            - "STANDARD": stores remarks that follow pdb standard
            - "NOT STANDARD": stores remarks not following the pdb standard
            - "ALL": stores all remarks
        "add_source_details" - Will add two remarks before the model tag: the path of the source file and 
            the original model number. 
        
    """
    keep_remarks = options.get_value("keep_remarks", default_value = "NONE")
    add_source_details = options.get_value("add_source_details", default_value = False)

    file_handler_out = open(out_pdb_name, "w")
    
    data = data_handler.get_data()
    
    merged_structure = data.get_all_elements()
    
    file_handler_out.write("REMARK 000 File created using Prody and pyProCT\n")
    
    if not keep_remarks == "NONE" and not add_source_details:
        prody.writePDBStream(file_handler_out, merged_structure, csets =  elements)
    else:
        all_remarks = filter_remarks(data.get_all_remarks(), subset= keep_remarks)
        all_model_numbers = data.get_all_model_numbers()
        
        current_model = 0
        for i, element_id in enumerate(elements): 
            if keep_remarks:
                remarks = all_remarks[element_id]
                file_handler_out.write("".join(remarks))
            
            if add_source_details:
                model_number = all_model_numbers[element_id]
                conf_source = data_handler.get_source_of_element(element_id).get_path()
                file_handler_out.write("REMARK source            : %s\n"%conf_source)
                file_handler_out.write("REMARK original model nr : %d\n"%model_number)
                file_handler_out.write("REMARK cluster id : %s\n"%ids[i])
                file_handler_out.write("REMARK cluster element : %d\n"%element_id)
            
            file_handler_out.write("MODEL"+str(current_model).rjust(9)+"\n")
            pdb_handler = cStringIO.StringIO()
            prody.writePDBStream(pdb_handler, merged_structure, csets=  element_id)
            # skip the first remark if any
            lines = filter(lambda line: line[0:6]!="REMARK" and line[0:5]!="MODEL" and line[0:6]!="ENDMDL", 
                           pdb_handler.getvalue().splitlines(True))
            pdb_handler.close()
            file_handler_out.write("".join(lines))
            file_handler_out.write("ENDMDL\n")
            current_model+=1

    file_handler_out.close()


