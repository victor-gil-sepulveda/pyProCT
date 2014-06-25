'''
Created on 25/06/2014

@author: victor
'''



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


