
from pyproclust.matrix.matrixCommon import generate_row
from pyproclust.matrix.completeMatrix import CompleteDistanceMatrix

class XPMConverter:
    """
    Reads an xpm file from g_rms and creates a complete matrix representation.
    Script de: 
        Tsjerk Wassenaar (tsjerkw at gmail.com) 
    Modificaciones y transformacion en clase de: 
        Victor Gil Sepulveda  (victor.gil.sepulveda at gmail.com)
    """
    def XPMConverter(self):
        pass
    
    def unquote(self,s):
        return s[1+s.find('"'):s.rfind('"')]
    
    def uncomment(self,s):
        return s[2+s.find('/*'):s.rfind('*/')]
    
    def col(self,c):
        color = c.split('/*')
        value = self.unquote(color[1])
        color = self.unquote(color[0]).split()
        return color[0], value
    
    def convert(self, xpm_file_handler_in):
        
        # Read in lines until we find the start of the array
        meta = [xpm_file_handler_in.readline()]
        while not meta[-1].startswith("static char *gromacs_xpm[]"):
            meta.append(xpm_file_handler_in.readline())
        
        # The next line will contain the dimensions of the array
        dim = xpm_file_handler_in.readline()
        
        # There are four integers surrounded by quotes
        nx, ny, nc, nb = [int(i) for i in self.unquote(dim).split()] #@UnusedVariable
        
        # The next dim[2] lines contain the color definitions
        # Each pixel is encoded by dim[3] bytes, and a comment
        # at the end of the line contains the corresponding value
        colors = dict([self.col(xpm_file_handler_in.readline()) for i in range(nc)])
        
        matrix_data = []
        for i in xpm_file_handler_in:
            if i.startswith("/*"):
                continue
            j = self.unquote(i)
            z = [colors[j[k:k+nb]] for k in range(0,nx,nb)]
            row_string = " ".join(z)+"\n"
            row = []
            generate_row(row_string, row)
            matrix_data.append(row)
        xpm_file_handler_in.close()
    
        return CompleteDistanceMatrix(matrix_data)

