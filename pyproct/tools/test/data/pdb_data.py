"""
Created on 28/02/2012

@author: victor
"""
pdb_1_num_of_models = 7
pdb1_num_of_atoms = 5
pdb_1_file_content = """MODEL        0
ATOM      3  CA  ILE     3      -0.039   0.638   3.156  1.00  1.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  1.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  1.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  1.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00
ENDMDL
MODEL        1
ATOM      3  CA  ILE     3      -1.039   0.638   3.156  1.00  2.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  2.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  2.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  2.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -2.039   0.638   3.156  1.00  3.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  3.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  3.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  3.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  3.00
ENDMDL
MODEL        3
ATOM      3  CA  ILE     3      -3.039   0.638   3.156  1.00  4.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  4.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  4.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  4.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  4.00
ENDMDL
MODEL        4
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  5.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  5.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  5.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  5.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  5.00
ENDMDL
MODEL        5
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  6.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  6.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  6.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  6.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  6.00
ENDMDL
MODEL        6
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  7.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  7.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  7.00
ENDMDL
"""

pdb_1_sub2_file_content = """MODEL 0
ATOM      3  CA  ILE     3      -0.039   0.638   3.156  1.00  1.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  1.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  1.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  1.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00
TER
MODEL 1
ATOM      3  CA  ILE     3      -1.039   0.638   3.156  1.00  2.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  2.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  2.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  2.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00
TER
MODEL 2
ATOM      3  CA  ILE     3      -2.039   0.638   3.156  1.00  3.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  3.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  3.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  3.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  3.00
TER
MODEL 3
ATOM      3  CA  ILE     3      -3.039   0.638   3.156  1.00  4.00
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
REMARK 400 REGULAR REMARK
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
ATOM      3  CA  ILE     3      -2.039   0.638   3.156  1.00  3.00
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
ATOM      3  CA  ILE     3      -0.039   0.638   3.156  1.00  1.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  1.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  1.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  1.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -2.039   0.638   3.156  1.00  3.00
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
REMARK 400 REGULAR REMARK
LOL LINE ------
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  4.00
ATOM      4  CA  PHE     4      -8.605   2.189  -0.292  1.00  4.00
ATOM      5  CA  VAL     5      -5.540   2.776  -2.541  1.00  4.00
ATOM      6  CA  LYP     6      -5.228   5.331  -5.451  1.00  4.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  4.00
ENDMDL
"""

premerged_pdb_1 = """MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  1.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  2.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00
ENDMDL
"""

premerged_pdb_2 = """MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  3.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  3.00
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  4.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  4.00
ENDMDL
"""

premerged_pdb_3 = """MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  5.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  5.00
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  6.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  6.00
ENDMDL
"""

merged_pdb = """MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  1.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  2.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00
ENDMDL
MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  3.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  3.00
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  4.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  4.00
ENDMDL
MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  5.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  5.00
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  6.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  6.00
ENDMDL
"""

merged_renumbered_pdb = """MODEL        0
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  1.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00
ENDMDL
MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  2.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  3.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  3.00
ENDMDL
MODEL        3
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  4.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  4.00
ENDMDL
MODEL        4
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  5.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  5.00
ENDMDL
MODEL        5
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  6.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  6.00
ENDMDL
"""

merged_1_5 = """MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  2.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00
ENDMDL
MODEL        5
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  6.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  6.00
ENDMDL
"""

merged_1_5_correlative = """MODEL        1
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  2.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  6.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  6.00
ENDMDL
"""

proto_pdb = """MODEL      100
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  1.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  1.00
ENDMDL
MODEL       26
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  2.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  2.00
ENDMDL
MODEL       48
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  3.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  3.00
ENDMDL
MODEL        2
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  4.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  4.00
ENDMDL
MODEL       56
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  5.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  5.00
ENDMDL
"""

proto_48_pdb = """MODEL       48
ATOM      3  CA  ILE     3      -9.039   0.638   3.156  1.00  3.00
ATOM      7  CA  THR     7      -3.037   3.652  -8.112  1.00  3.00
ENDMDL
"""

chain_padding_proto_1 = """MODEL        1
ATOM      1  CA  ILE A   1       1.000   2.000   3.000  0.00  0.00
ATOM      2  CA  ILE B   2       4.000   5.000   6.000  0.00  0.00
ATOM      3  CA  ILE B   2       7.000   8.000   9.000  0.00  0.00
ATOM      4  CA  ILE C   3      10.000  11.000  12.000  0.00  0.00
ATOM      5  CA  ILE C   3      13.000  14.000  15.000  0.00  0.00
ATOM      6  CA  ILE C   3      16.000  17.000  18.000  0.00  0.00
ENDMDL
"""

chain_padding_proto_2 = """MODEL        1
ATOM      1  CA  ILE A   1       1.000   2.000   3.000  0.00  0.00
ATOM      2  CA  ILE C   2       4.000   5.000   6.000  0.00  0.00
ATOM      3  CA  ILE C   2       7.000   8.000   9.000  0.00  0.00
ATOM      4  CA  ILE C   2      10.000  11.000  12.000  0.00  0.00
ATOM      5  CA  ILE B   3      13.000  14.000  15.000  0.00  0.00
ATOM      6  CA  ILE B   3      16.000  17.000  18.000  0.00  0.00
ENDMDL
MODEL        2
ATOM      1  CA  ILE A   1      19.000  20.000  21.000  0.00  0.00
ATOM      2  CA  ILE C   2      22.000  23.000  24.000  0.00  0.00
ATOM      3  CA  ILE C   2      25.000  26.000  27.000  0.00  0.00
ATOM      4  CA  ILE C   2      28.000  29.000  30.000  0.00  0.00
ATOM      5  CA  ILE B   3      31.000  32.000  33.000  0.00  0.00
ATOM      6  CA  ILE B   3      34.000  35.000  36.000  0.00  0.00
ENDMDL
"""

switched_pdb_data = """MODEL        1
ATOM      1  ATO FAK A  00       1.000   2.000   3.000  0.00  0.00
ATOM      1  Pt1 CPT X  25     -33.115   1.294  -1.163  0.00  0.00
ENDMDL
MODEL        2
ATOM      1  Pt1 CPT X  25     -32.555  -2.500  -5.367  0.00  0.00
ATOM      1  ATO FAK A  00       4.000   5.000   6.000  0.00  0.00
ENDMDL
MODEL        3
ATOM      1  ATO FAK A  00       7.000   8.000   9.000  0.00  0.00
ATOM      1  Pt1 CPT X  25     -33.257   5.280  -8.441  0.00  0.00
ENDMDL
MODEL        4
ATOM      1  Pt1 CPT X  25     -32.306   6.517  -1.544  0.00  0.00
ATOM      1  ATO FAK A  00      10.000  11.000  12.000  0.00  0.00
ENDMDL
MODEL        5
ATOM      1  ATO FAK A  00      13.000  14.000  15.000  0.00  0.00
ATOM      1  Pt1 CPT X  25     -30.494  10.390  -3.066  0.00  0.00
ENDMDL
"""