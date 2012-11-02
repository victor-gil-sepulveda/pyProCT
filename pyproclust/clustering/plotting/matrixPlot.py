'''
Created on 31/01/2012

@author: victor
'''
import numpy as np
import colorsys as color
import matplotlib
import matplotlib.pyplot
from pyproclust.matrix.completeMatrix import CompleteDistanceMatrix

def hvs_to_rgb_converter_yellow_to_red(norm_value):
    """
    Example of 'color descriptor' to be used with the matrix to image functions.
    It has to return the color value associated to a normalized value of the distance
    matrix. This one establishes a relationship between a [0,1] value and a RGB value
    using a yellow to red palette.
    """
    point_color = np.array(color.hsv_to_rgb(norm_value*45/360,0.8,0.8))*255
    return (int(point_color[0]),int(point_color[1]),int(point_color[2]))


def __putpixel(norm_value, pixel_col, i,j,color_descriptor):
    """
    Puts a pixel into a pixel collection.
    - Getting the pixel collection:
    from PIL import Image
    image = Image.new('RGB', (row_len,row_len))
    pix = image.load()
    pix[j,i] = ...
    
    """
    #point_color = np.array(color.hsv_to_rgb(norm_value*45/360,0.8,0.8))*255
    pixel_col[i,j] = color_descriptor(norm_value)



def condensed_matrix_to_image(condensed_matrix, image_pixels, color_descriptor=hvs_to_rgb_converter_yellow_to_red):
    """
    Writes the contents of a NORMALIZED[0,1] condensed distance matrix in a pixel collection which 
    offers[x,y] access (as those in PIL).  
    """
    #image = Image.new('RGB', (row_len,row_len))
    #pix = image.load()
    row_lenght = condensed_matrix.row_length
    
    for i in range(0,row_lenght):
        for j in range(0,row_lenght):
            norm_value = condensed_matrix[i,j]
            # The image is transposed
            __putpixel(norm_value, image_pixels, j, i,color_descriptor)


def data_matrix_to_image(matrix,image_pixels, color_descriptor):
    """
    Writes the contents of a NORMALIZED[0,1] distance matrix in a pixel collection which 
    offers[x,y] access (as those in PIL).  
    """
    column_len = len(matrix)
    row_len = len(matrix[0])

    #image = Image.new('RGB', (row_len,column_len))
    #pix = image.load()
    for i in range(column_len):
        for j in range(row_len):
            # The image is transposed
            #point_color = np.array(color.hsv_to_rgb(matrix[i][j]*45/360,0.8,0.8))*255
            #pix[j,i] = (int(point_color[0]),int(point_color[1]),int(point_color[2]))
            
            __putpixel(matrix[i][j],image_pixels,j,i,color_descriptor)
    #image.save(file_name) 


def plot_matrix( matrix, save_path = None):
    """
    This function plots a SQUARE distance matrix using matplotlib and writes it to disk
    if a save_path is given.
    """
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(matrix.get_data(),norm = matplotlib.colors.Normalize())
    im.set_cmap("hot")
    cbar = fig.colorbar(im, ticks=[-1, 0, 1])
    cbar.ax.set_yticklabels(['< -1', '0', '> 1'])
    matplotlib.pyplot.show()
    if save_path:
        fig.savefig("lol.png",format="png")
        
if __name__ == '__main__':
    matrix = [[100, 100, 67.37, 23.53],
              [100, 100, 92.68, 31.02],
              [67.37, 92.68, 100, 23.53],
              [23.53, 31.02, 23.53, 100]]
    matrix = CompleteDistanceMatrix(matrix)
    plot_matrix(matrix,"graph")