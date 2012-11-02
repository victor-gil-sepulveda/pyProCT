'''
Created on 28/02/2012

@author: victor
'''
pdb_1_num_of_models = 7
pdb_1_file_content = """MODEL 0
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  1.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  1.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  1.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  1.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00            
TER           
MODEL 1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  2.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  2.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  2.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  2.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00            
TER            
MODEL 2        
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  3.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  3.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  3.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  3.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  3.00            
TER            
MODEL 3
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  4.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  4.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  4.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  4.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  4.00            
TER
MODEL 4
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  5.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  5.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  5.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  5.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  5.00            
TER            
MODEL 5
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  6.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  6.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  6.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  6.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  6.00            
TER
MODEL 6
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  7.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  7.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  7.00            
TER     
"""

pdb_2_num_of_models = 4
pdb_2_file_content = """MODEL 0
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  1.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  1.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  1.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  1.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00            
ENDMDL
THINGHIE ************
REMARK 0 this 
REMARKS 0 is
REMARKS 0 a remark 
MODEL 1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  2.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  2.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  2.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  2.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00            
ENDMDL   
LOL LINE ------
REMARKS 0 this 
REMARKS 0 is
REMARKS 0 a remark         
MODEL 2        
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  3.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  3.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  3.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  3.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  3.00            
ENDMDL
THINGHIE ************
THINGIE *******
LOL LINE ------            
MODEL 3
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  4.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  4.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  4.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  4.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  4.00            
ENDMDL
"""

extracted_pdbs_1 = """MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  3.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  3.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  3.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  3.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  3.00            
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  6.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  6.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  6.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  6.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  6.00            
ENDMDL
"""

extracted_pdbs_2 = """MODEL        0
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  1.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  1.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  1.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  1.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00            
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  3.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  3.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  3.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  3.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  3.00            
ENDMDL
"""

extracted_pdbs_3 = """MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  2.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  2.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  2.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  2.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00            
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  4.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  4.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  4.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  4.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  4.00            
ENDMDL
"""

extracted_pdbs_4 = """THINGHIE ************
REMARK 0 this 
REMARKS 0 is
REMARKS 0 a remark 
MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  2.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  2.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  2.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  2.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00            
ENDMDL
THINGHIE ************
THINGIE *******
LOL LINE ------            
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  4.00            
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  4.00            
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  4.00            
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  4.00            
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  4.00            
ENDMDL
"""