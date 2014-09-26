"""
Created on 10/02/2014

@author: victor
"""
import os.path
import prody

class RepresentativesPostAction(object):

    KEYWORD = "representatives"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, data_handler, workspaceHandler, matrixHandler, generatedFiles):

        medoids = clustering.get_medoids(matrixHandler.distance_matrix)
    
        pdb_name = postprocessing_parameters.get_value("filename", default_value = "representatives")
        
        representatives_file_path = os.path.join( workspaceHandler["results"],"%s.pdb"%pdb_name)
    
        save_representatives( medoids,
                              representatives_file_path,
                              data_handler,
                              postprocessing_parameters)
    
        generatedFiles.append({
                               "description":"Cluster representatives",
                               "path":os.path.abspath(representatives_file_path),
                               "type":"pdb"
        })

def save_representatives(representatives,
                         out_pdb_name,
                         data_handler,
                         options):
    """
    Saves a pdb file containing the most representative elements of the clustering.

    @param representatives: A list of the representative elements of the clustering we want to extract.

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
    
    if not keep_remarks == "NONE" and not add_source_details:
        prody.writePDBStream(file_handler_out, merged_structure, csets =  representatives)
    else:
        all_remarks = data.get_all_remarks()
        all_model_numbers = data.get_all_model_numbers()
        
        for element_id in representatives: 
            if keep_remarks:
                remarks = all_remarks[element_id]
                file_handler_out.write("".join(remarks))
            
            if add_source_details:
                model_number = all_model_numbers[element_id]
                conf_source = data_handler.get_source_of_element(element_id).get_path()
                file_handler_out.write("REMARK source            : %s\n"%conf_source)
                file_handler_out.write("REMARK original model nr : %d\n"%model_number)
            
            prody.writePDBStream(file_handler_out, merged_structure, csets =  element_id)
    
    file_handler_out.close()


