"""
Created on 23/02/2012

@author: victor
"""
import random
import numpy
from PIL import Image, ImageDraw
import scipy.spatial.distance as distance
from pyRMSD.condensedMatrix import CondensedMatrix

def create_matrices(data, verbose = False):
    """
    Creates all the matrices for the observations (datasets) and returns both.
    """
    condensed_matrices = {}
    all_observations = {}
    for dataset_name in data.all_datasets:
        dataset = data.all_datasets[dataset_name]
        # Creating the matrix
        observations = dataset_loading_2D(dataset,data.scale_factor[dataset_name])

        all_observations[dataset_name] = observations
        condensed_matrix_data = distance.pdist(observations)
        condensed_matrix = CondensedMatrix(condensed_matrix_data)
        condensed_matrices[dataset_name] = condensed_matrix
        if verbose:
            print "Matrix for %s:"%dataset_name
            print "-----------------------"
            print "Max dist. = ",condensed_matrix.calculateMax()
            print "Min dist. = ",condensed_matrix.calculateMin()
            print "Mean dist. = ",condensed_matrix.calculateMean()
            print "Variance = ",condensed_matrix.calculateVariance()
            print "-----------------------\n"
    return condensed_matrices, all_observations

def params_to_string(params):
    """
    Converts a params dictionary in a string suitable for file name creation.

    @param params: The clustering algorithm params (kwargs) to convert.

    @return: Stringification of this params.
    """
    s = ""
    for k in params:
        try:
            value = params[k]
            if isinstance( value, ( int, long ) ):
                s += "%s_%d_"%(k, value)
            else:
                s += "%s_%.3f_"%(k, value)
        except TypeError:
            s += "%s_%s_"%(k, params[k])
    return s[:-1]

def dataset_loading_2D(dataset_string, scale_factor = 1):
    """
    Creates a list of observations from a 2D dataset.
    """
    observations = []
    lines = dataset_string.split("\n")
    for l in lines:
        try:
            nums = l.split()
            observations.append((float(nums[0])*scale_factor,float(nums[1])*scale_factor))
        except IndexError:
            pass
        except:
            print "[Error dataset_loading_2D] Impossible to parse line. Exiting..."
            exit()
    return observations

def get_2D_bounding_box(observations):
    """
    Gets the maximum and minimum x and y in order to construct the bounding box.
    """
    x = []
    y = []
    # slow but elegant! :P
    for o in observations:
        x.append(o[0])
        y.append(o[1])
    return (min(x),min(y)),(max(x),max(y))

def create_canvas(dataset_observations, scale, margin):
    """
    Creates the canvas
    """
    boundaries = get_2D_bounding_box(dataset_observations)
    width = boundaries[1][0]-boundaries[0][0]
    height =  boundaries[1][1]-boundaries[0][1]
    image = Image.new('RGB', (int(width*scale+2*margin),int(height*scale+2*margin)))
    draw = ImageDraw.Draw(image)

    return (image,draw,boundaries)

def draw_point_into_canvas(dataset_observations, i, draw, scale, margin, color, boundaries, print_number, shape_function = None):
    """
    Draws dataset observation 'i' in the canvas, with color = 'color'
    """
    position = ((dataset_observations[i][0] - boundaries[0][0])*scale+margin,(dataset_observations[i][1] - boundaries[0][1])*scale+margin)

    if print_number:
        draw.text(position, str(i),fill=color)

    if shape_function is None :
        draw.point(position)
    else:
        boundaries = (position[0]-int((0.5*scale)),position[1]-int((0.5*scale)),
                    position[0]+int((0.5*scale)),position[1]+int((0.5*scale)))
        shape_function(draw, boundaries, color)
    return position

def show_2D_dataset(dataset_observations, scale, margin = 0, cutoff = 0, cutoffed_nodes = []):
    """
    Draws the dataset into a PIL image and returns it.
    The scale parameter is a pixel multiplication factor for the image size, the margin
    is a quantity in pixels that will be added to each of the sides of the bounding
    box (previous to the scale multiplication). The cutoff is the radius used to do the gromos
    clusterization, and the cutoffed_nodes list stores the nodes we want to draw their cutoffs.
    """
    scaled_radi = cutoff*scale

    (image,draw,boundaries) = create_canvas(dataset_observations, scale, margin)

    for i in range(len(dataset_observations)):
        position = draw_point_into_canvas(dataset_observations,i,draw,scale,margin,"#ff0000",boundaries)
        if i in cutoffed_nodes:
            draw.ellipse((position[0]-scaled_radi,position[1]-scaled_radi,\
                          position[0]+scaled_radi,position[1]+scaled_radi))
    return image

def generate_color_list(number_of_colors):
    """
    Creates a number of random colors
    """
    color_list = []
    for i in range(number_of_colors): #@UnusedVariable
        r = random.randint(40, 255)
        g = random.randint(40, 255)
        b = random.randint(40, 255)
        scolor = ("#%02x%02x%02x"%(r,g,b))
        color_list.append(scolor)
    return color_list

def generate_faded_red(alpha):
    """
    Creates a number of random colors
    """
    if alpha > 0 and alpha < 0.2:
        alpha = 0.2 # helps to visualize weak relationships

    r = int(alpha * 255.)
    g = 0.0
    b = 0.0
    return  "#%02x%02x%02x"%(r,g,b)

def draw_cross(draw, boundaries, color):
    xsup = boundaries[0]
    ysup = boundaries[1]
    xinf = boundaries[2]
    yinf = boundaries[3]
    draw.line([(xsup, ysup), (xinf, yinf)], width = 2, fill = color)
    draw.line([(xsup, yinf), (xinf, ysup)], width = 2, fill = color)

def draw_circle(draw, boundaries, color):
    draw.ellipse(boundaries, fill = color)

def draw_square(draw, boundaries, color):
    draw.rectangle(boundaries, fill = color)

def draw_rombo(draw, boundaries, color):
    xsup = boundaries[0]
    ysup = boundaries[1]
    xinf = boundaries[2]
    yinf = boundaries[3]

    draw.polygon([(xsup, ysup+0.5*(yinf-ysup)),
                  (xsup+0.5*(xinf-xsup), ysup),
                  (xinf, ysup+0.5*(yinf-ysup)),
                  (xsup+0.5*(xinf-xsup),yinf)], fill = color)

def draw_triangle(draw, boundaries, color):
    xsup = boundaries[0]
    ysup = boundaries[1]
    xinf = boundaries[2]
    yinf = boundaries[3]
    draw.polygon([(xsup,yinf),(xsup+0.5*(xinf-xsup),ysup),(xinf,yinf)], fill = color)

def show_2D_dataset_clusters(dataset_observations, clusterization, scale, margin = 0, print_numbers = False):
    """
    Generates an image with a 2D dataset drawn, where alll the points belonging to the same
    cluster have the same color, different from the color of the other clusters.
    """
    (image,draw,boundaries) = create_canvas(dataset_observations, scale, margin)

    color_list = generate_color_list(len(clusterization.clusters))
    available_shape_functions = [draw_cross, draw_circle, draw_square, draw_triangle, draw_rombo]
    for i, c in enumerate(clusterization.clusters):
        color = color_list.pop()
        for j in c.all_elements:
            draw_point_into_canvas(dataset_observations, j, draw, scale, margin, color, boundaries, print_numbers, available_shape_functions[i%5])
    return image

def generate_similarity_network(W, dataset_observations, scale, margin = 0, print_numbers = False):
    norm_W = W / numpy.max(W)

    (image, draw, boundaries) = create_canvas(dataset_observations, scale, margin)

    for i in range(len(dataset_observations)-1):
        i_position = ((dataset_observations[i][0] - boundaries[0][0])*scale+margin,(dataset_observations[i][1] - boundaries[0][1])*scale+margin)
        for j in range(i+1, len(dataset_observations)):
            if norm_W[i][j] > 0.0:
                j_position = ((dataset_observations[j][0] - boundaries[0][0])*scale+margin,(dataset_observations[j][1] - boundaries[0][1])*scale+margin)
                #print generate_faded_red(norm_W[i][j])
                draw.line([i_position, j_position], width = 2, fill = generate_faded_red(norm_W[i][j]))

    for i in range(len(dataset_observations)):
        position = ((dataset_observations[i][0] - boundaries[0][0])*scale+margin,(dataset_observations[i][1] - boundaries[0][1])*scale+margin)
        draw.point(position, fill="#FFFFFF")
        if print_numbers:
            draw.text(position, str(i), fill="#FFFFFF")

    return image
