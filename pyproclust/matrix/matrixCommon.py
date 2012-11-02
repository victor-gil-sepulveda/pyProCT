'''
Created on 12/03/2012

@author: victor
'''

def generate_row(line,row):
    """
    Creates a row from a line of text.
    """    
    numbers = line.split()
    for n in numbers:
        row.append(float(n))
        
def assertMatrixAlmostEqual(test,matrix1_data,matrix2_data,tol):
    """
    Tests if a matrix is equal to another within a given precission (to be used with
    the unit testing framework).
    """
    if(len(matrix1_data)!= len(matrix2_data)):
        return False
    for i in range(len(matrix1_data)):
        for j in range(len(matrix1_data[0])):
            test.assertAlmostEqual(matrix1_data[i][j], matrix2_data[i][j], tol)
    return True