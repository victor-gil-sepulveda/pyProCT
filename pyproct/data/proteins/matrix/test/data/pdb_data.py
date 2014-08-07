"""
Created on 25/06/2014

@author: victor
"""



pdb_data = """MODEL        1
ATOM      1  ATO FAK A  00       1.000   2.000   3.000  0.00  0.00
ATOM      1  Pt1 CPT X  25     -33.115   1.294  -1.163  0.00  0.00
ENDMDL
MODEL        2
ATOM      1  ATO FAK A  00       4.000   5.000   6.000  0.00  0.00
ATOM      1  Pt1 CPT X  25     -32.555  -2.500  -5.367  0.00  0.00
ENDMDL
MODEL        3
ATOM      1  ATO FAK A  00       7.000   8.000   9.000  0.00  0.00
ATOM      1  Pt1 CPT X  25     -33.257   5.280  -8.441  0.00  0.00
ENDMDL
MODEL        4
ATOM      1  ATO FAK A  00      10.000  11.000  12.000  0.00  0.00
ATOM      1  Pt1 CPT X  25      32.306   6.517  -1.544  0.00  0.00
ENDMDL
MODEL        5
ATOM      1  ATO FAK A  00      13.000  14.000  15.000  0.00  0.00
ATOM      1  Pt1 CPT X  25      30.494  10.390  -3.066  0.00  0.00
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

chain_padding_proto_3 = """MODEL        1
ATOM      1  CA  ILE A   1       1.000   2.000   3.000  0.00  0.00
ATOM      2  CA  ILE B   2       4.000   5.000   6.000  0.00  0.00
ATOM      3  CA  ILE B   2       7.000   8.000   9.000  0.00  0.00
ATOM      4  CA  ILE C   3      10.000  11.000  12.000  0.00  0.00
ATOM      5  CA  ILE C   3      13.000  14.000  15.000  0.00  0.00
ATOM      6  CA  ILE C   3      16.000  17.000  18.000  0.00  0.00
ATOM      7  CA  ILE D   4      19.000  20.000  21.000  0.00  0.00
ATOM      8  CA  ILE D   4      22.000  23.000  24.000  0.00  0.00
ATOM      9  CA  ILE E   5      25.000  26.000  27.000  0.00  0.00
ATOM     10  CA  ILE E   5      28.000  29.000  30.000  0.00  0.00
ATOM     11  CA  ILE E   5      31.000  32.000  33.000  0.00  0.00
ENDMDL
MODEL        2
ATOM      1  CA  ILE A   1      41.000  42.000  43.000  0.00  0.00
ATOM      2  CA  ILE B   2      44.000  45.000  46.000  0.00  0.00
ATOM      3  CA  ILE B   2      47.000  48.000  49.000  0.00  0.00
ATOM      4  CA  ILE C   3      50.000  51.000  52.000  0.00  0.00
ATOM      5  CA  ILE C   3      53.000  54.000  55.000  0.00  0.00
ATOM      6  CA  ILE C   3      56.000  57.000  58.000  0.00  0.00
ATOM      7  CA  ILE D   4      59.000  60.000  61.000  0.00  0.00
ATOM      8  CA  ILE D   4      62.000  63.000  64.000  0.00  0.00
ATOM      9  CA  ILE E   5      65.000  66.000  67.000  0.00  0.00
ATOM     10  CA  ILE E   5      68.000  69.000  70.000  0.00  0.00
ATOM     11  CA  ILE E   5      71.000  72.000  73.000  0.00  0.00
ENDMDL
"""


