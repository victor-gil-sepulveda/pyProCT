"""
Created on 09/12/2013

@author: victor
"""

if __name__ == '__main__':
    output = open("fake.pdb","w")

    for i in range(1,1000):
        output.write("MODEL"+str(i).rjust(9)+"\nATOM      1  ATO FAK A  00       1.000   2.000   3.000  0.00  0.00\nENDMDL\n")
    output.close()