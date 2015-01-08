"""
Created on 15/10/2012

@author: victor
"""
import numpy
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw
import matplotlib.pyplot as plt
import pylab
import PIL.ImageFont as ImageFont
from pyRMSD.condensedMatrix import CondensedMatrix
import math

def shrink_matrix(this_matrix, to_have_this_size):
    """
    Uses a stencil technique to scale down the matrix in order to be able to represent it in an image file
    [Warning] It does not any size check. If the matrix is already smaller than the expected size, the
    function behaviour is undefined.
    
    :param this_matrix: The matrix we want to reduce.
    :param to_have_this_size: The new side length of the matrix.
    
    :return: A resized matrix.
    """
    max_dim = to_have_this_size
    row_dim = this_matrix.row_length

    box_size = int(math.ceil(row_dim/max_dim))
    number_of_boxes = int(math.ceil(row_dim /box_size))
    box_mid = int(math.ceil(box_size/2.))
    max_dim = number_of_boxes * box_size
    
    print "* Producing reduced matrix image (from %d pixels to %d pixels)."%(row_dim, to_have_this_size)
    print "* Box_size: %d pixels"%box_size

    tmp_condensed = CondensedMatrix(numpy.zeros(int((number_of_boxes*(number_of_boxes-1))/2.), dtype=numpy.float))

    for k in range(0, max_dim, box_size):
        for l in range(0, max_dim, box_size):
            values = []
            for i in range(-box_mid, box_mid):
                for j in range(-box_mid, box_mid):
                    if(k+i > 0 and l+j > 0 and k+i < max_dim and l+j < max_dim):
                        if k+i == l+j:
                            values.append(0)
                        else:
                            values.append(this_matrix[k+i,l+j])
            if(k/box_size != l/box_size):
                if values != []:
                    tmp_condensed[int(k/box_size), int(l/box_size)] = numpy.mean(values)
                else:
                    tmp_condensed[int(k/box_size), int(l/box_size)] = 0
    return tmp_condensed

def matrixToImage(condensed_distance_matrix, matrix_image_file, max_dim = 1000, diagonal_value = 0, observer = None):
    """
    Generates a plot of the distance matrix given as argument and stores it into disk.

    @param condensed_distance_matrix: Is the matrix (CondensedMatrix) from which we want to create
    the image.

    @param matrix_image_file: The path of the image file to create (with file extension).

    @param max_dim: Maximum dimensions of the image. If the matrix is bigger it is rescaled.
    """

    if condensed_distance_matrix.row_length > max_dim:
        if not observer is None:
            observer.notify("Matrix To Image","Reescale","Dimension of the matrix ("+str(condensed_distance_matrix.row_length)+") is bigger than"+
            " the maximum dimension for imaging ("+str(max_dim)+")")
        matrix = shrink_matrix(condensed_distance_matrix, max_dim)
    else:
        matrix = condensed_distance_matrix

    complete = numpy.zeros([ matrix.row_length]*2, dtype=numpy.float)

    # fill diagonal if needed
    if diagonal_value != 0:
        for i in range(matrix.row_length):
            complete[i][i] = diagonal_value

    # fill matrix
    for i in range(matrix.row_length-1):
        for j in range(i+1,matrix.row_length):
            complete[i][j] = matrix[i,j]
            complete[j][i] = matrix[i,j]

    plt.clf()
    if condensed_distance_matrix.row_length > max_dim:
        imgplot = plt.imshow(complete, interpolation='nearest')
    else:
        imgplot = plt.imshow(complete, interpolation='none')
    imgplot.set_cmap('gray')
    plt.colorbar()
    plt.savefig(matrix_image_file)
#     plt.show()
    plt.close()

def pieChartCreation(graph_size, fracs, name1, name2, colors):
    """
    Creates the big pie chart of the report. In this pie chart one fraction represents the amount of space
    sampled together by both trajectories, and the other two fractions represent the amount of space that was
    sampled by either A or B.

    @param graph_size: The graph size in pixels (?)
    @param fracs: a list or tuple containing the total number of elements of pure A, pure B, and mixed clusters.
    @param name1: String identifying the first trajectory.
    @param name2: String identifying the second trajectory.
    @param colors: dictionary with color descriptions for "A", "B" and "M" (Mixed) written in string format (for
    example "#FFFFFF" for white)

    @return : A PIL image of the pie chart.
    """
    all_labels = ["A","B","Mixed"]
    all_colors = [colors['A'], colors['B'],colors['M']]
    this_labels = []
    this_fracs = []
    this_colors = []
    mydpi = 100
    fig = pylab.figure(figsize=(int(graph_size[0]/mydpi),int(graph_size[1]/mydpi)),dpi = mydpi)
    fig.set_facecolor('white')
    for i in range(len(fracs)):
        if fracs[i]!=0:
            this_labels.append(all_labels[i])
            this_fracs.append(fracs[i])
            this_colors.append(all_colors[i])
    pylab.pie(this_fracs, labels=this_labels, autopct='%1.1f%%', shadow=False,colors=this_colors)
    pylab.title(shorten_name(name1)+" vs "+shorten_name(name2))
    return fig2img(fig)

def barGraphCreation(A_sizes, B_sizes, cluster_sizes, types, total_size, colors, graph_size):
    """
    Draws a bar graph with information about the number of elements in the first ten biggest clusters.

    @param A_sizes: Total number of elements in 'A' clusters, f.ex. [2,3,0,4,5] In this case the third cluster would be of type B.
    @param B_sizes: The same, for 'B' type clusters.
    @param cluster_sizes: List containing the size of each cluster.
    @param types: List containing the type of each cluster ('A','B' or 'M' (Mixed))
    @param total_size: Number of clustered elements.
    @param colors: dictionary with color descriptions for "A", "B" and "M" (Mixed) written in string format (for
    example "#FFFFFF" for white)
    @param graph_size: Size of the image (width, height)

    @return : A PIL image of the graph bar.
    """
    mydpi = 100.
    fig = plt.figure(figsize = (graph_size[0]/mydpi, (graph_size[1])/mydpi), dpi = mydpi)
    fig.set_facecolor('white')
    ax = fig.add_subplot(111)
    number_of_bars = min(10,len(A_sizes))
    A_values = numpy.array(A_sizes[0:number_of_bars])*100. /  total_size
    B_values = numpy.array(B_sizes[0:number_of_bars])* 100. / total_size

    # Add pure clusters
    for i in range(number_of_bars):
        if types[i]=='A':
            A_values[i] = cluster_sizes[i]*100./total_size
        if types[i]=='B':
            B_values[i] = cluster_sizes[i]*100./total_size

    x = numpy.array(range(len(A_values)))
    ax.bar(left = x,width = 1.0, height = A_values, color = colors['A'])
    ax.bar(left = x, width = 1.0, bottom = A_values, height = B_values, color = colors['B'])
    ax.set_ylabel('%')

    pylab.yticks(fontsize=8)
    pylab.xticks(fontsize=7)
    fig.text(0.75,0.8,"Acc: %0.2f"%(numpy.sum(cluster_sizes)*100/total_size),size= 7)

    return fig2img(fig)

def ballGraph(number_of_balls, size, max_box_diameter, h_ball_separation, v_ball_separation,\
                             sizes, alphas, colors, cluster_types, A_sizes, B_sizes):
    """
    It draws an image with N spheres (where N is typically the number of clusters we have). This spheres are
    outlined and in the case they're representing a 'M' type cluster, they have an inner circle representing the amount
    of elements of 'A' and 'B' that it has.

     @param number_of_balls: The number of data analysis we have done (usually number of clusters).
     @param size: Size of the image to hold the ball graph.
     @param max_box_diameter: Is the diameter of the biggest possible sphere.
     @param h_ball_separation: Separation of each sphere in the x axis.
     @param v_ball_separation: Separation of each sphere in the y axis (usually enough to fit the 'data card')
     @param sizes: List containing the size of each cluster.
     @param alphas: Is a list containing the dispersion for each cluster.
     @param colors: Dictionary of colors for each type ('A','B' and 'M') with RGB color in tuple form (f.ex. (255,255,255) is white)
     @param cluster_types: List containing the type of each cluster ('A','B' or 'M' (Mixed))
     @param A_sizes: Total number of elements in 'A' clusters, f.ex. [2,3,0,4,5] In this case the third cluster would be of type B.
     @param B_sizes: The same, for 'B' type clusters.

     @return: The PIL image with the drawing and the maximum calculated diameter of the biggest sphere.
    """

    # Refactorizable
    w,h = size
    max_number_of_balls_in_a_row  = 0
    max_number_of_balls_in_a_column = 0
    max_tmp_diameter = max_box_diameter - 4
    while  max_number_of_balls_in_a_row * max_number_of_balls_in_a_column <= number_of_balls:
        max_number_of_balls_in_a_row = int(w / (float(max_tmp_diameter)+h_ball_separation))
        max_number_of_balls_in_a_column = int(h / (float(max_tmp_diameter)+v_ball_separation))
        max_tmp_diameter -= 1

    # Normalize ball sizes and alphas (based on dispersion)
    norm_sizes = normalize(sizes,max_box_diameter)
    # Reverse values and normalize between 1 and 0.3 (0.3 is the alpha of the cluster with most dispersion)
    norm_alphas = normalize_in_range(1-numpy.array(normalize_in_range(alphas,0,1)),0.3,1)
    # Paint it!!!
    canvas = Image.new("RGBA", size,color=(255,)*4)
    k = 0
    for i in range(max_number_of_balls_in_a_column):
        for j in range(max_number_of_balls_in_a_row):
            k = i*max_number_of_balls_in_a_row+j
            if k < number_of_balls:
                this_diameter = norm_sizes[k]*float(max_tmp_diameter)/float(max_box_diameter)
                pos = ((j*(max_tmp_diameter+h_ball_separation)),\
                           (i*(max_tmp_diameter+v_ball_separation)))

                centering_offset = max_box_diameter/2. - this_diameter/2
                ball = drawBall((max_box_diameter,)*2, this_diameter, colors[cluster_types[k]],norm_alphas[k])
                if cluster_types[k] != 'M':
                    canvas.paste(ball,tuple_to_int(pos))
                else:
                    ball_with_inner_circles = circledPercentages(ball, (centering_offset,centering_offset), this_diameter, 0.3, A_sizes[k], B_sizes[k], colors)
                    canvas.paste(ball_with_inner_circles,tuple_to_int(pos))
    return canvas, max_tmp_diameter

def drawBall(canvas_size, ball_diameter, fill_color, alpha):
    """
    Draws the image of black-outlined sphere representing a cluster.

    @param ball_diameter: The diameter of the sphere...
    @param fill_color: Color used to fill the sphere.
    @param alpha: Is the alpha value of the image (usually the

    @return: A PIL image of the resulting sphere.
    """
    canvas = Image.new("RGBA", canvas_size, color=(255,)*4)
    draw = ImageDraw.Draw(canvas)
    centering_offset = int(max(canvas_size)/2. - ball_diameter/2)
    position = (centering_offset,centering_offset,centering_offset+ball_diameter,centering_offset+ball_diameter)
    draw.ellipse(position,fill=fill_color)
    blending_canvas = Image.new("RGBA", canvas_size,color=(255,)*4)
    border_canvas = Image.new("RGBA", canvas_size,color=(255,)*4)
    drawCircleBorder(border_canvas, position, 2)
    return Image.composite(Image.blend(blending_canvas,canvas,alpha), border_canvas, border_canvas)

def drawCircleBorder(canvas, bbox, width):
    """
    Draws a black border over a sphere PIL image.

    @param canvas: The PIL image with the sphere.
    @param bbox: The description of this image's bounding box in coordinates (upper left x, ul y, bottom right x, br y).
    @param width: The width of the border.
    """
    draw = ImageDraw.Draw(canvas)
    border_position = (bbox[0]-width,bbox[1]-width,bbox[2]+width,bbox[3]+width)
    draw.ellipse(border_position,fill=(0,)*4)
    draw.ellipse(tuple_to_int(bbox),fill=(255,)*4)

def circledPercentages(ball_canvas, position, ball_diameter, radius_percent, sizeA, sizeB, colors):
    """
    Draws the inner part of the sphere (if it's of type 'M') with the percentage of 'A' and 'B' parts.

    @param ball_canvas: Sphere image to draw the inner part.
    @param position: The corner of the bounding box of the sphere, relative to the image.
    @param ball_diameter: The diameter of the sphere drawn in ball_canvas.
    @param radius_percent: Is the percent of the radius dedicated to the percentage representation.
    @param sizeA: Number of elements of 'A' in the 'M' cluster.
    @param sizeB: Number of elements of 'B' in the 'M' cluster.
    @param colors: Dictionary of colors for each type ('A','B' and 'M') with RGB color in tuple form (f.ex. (255,255,255) is white)

    @return: A PIL image with the complete sphere drawing (the complete representation of an 'M' cluster).
    """
    ball_radius = ball_diameter/2.
    width = radius_percent * ball_radius

    inner_bbox = (position[0]+0.5*ball_radius, position[1]+0.5*ball_radius,\
                  position[0]+ball_diameter-(0.5*ball_radius), position[1]+ball_diameter-(0.5*ball_radius))

    outer_bbox = (inner_bbox[0]-width, inner_bbox[1]-width,\
                  inner_bbox[2]+width, inner_bbox[3]+width)

    A_mask = Image.new("RGBA", ball_canvas.size,color=(255,)*4)
    drawCircleBorder(A_mask,tuple_to_int(inner_bbox),width)
    A_color = Image.new("RGBA", ball_canvas.size,color=colors['A'])
    A_circled =  Image.composite(ball_canvas, A_color, A_mask)

    B_mask = Image.new("RGBA", ball_canvas.size,color=(255,)*4)
    drawCircleBorder(B_mask,tuple_to_int(inner_bbox),width)
    B_color = Image.new("RGBA", ball_canvas.size,color=colors['B'])
    draw = ImageDraw.Draw(B_mask)
    B_percent = sizeB / float(sizeA+sizeB)
    draw.pieslice(tuple_to_int(outer_bbox), 0, int(360 *(1- B_percent)), fill = (255,)*4)

    # Separation line from A-B
    draw = ImageDraw.Draw(B_color)
    draw.pieslice(tuple_to_int(outer_bbox), -1, int(360 *(1- B_percent))+1, outline = (0,)*4 )

    AB_circled = Image.composite(A_circled, B_color, B_mask)
    outer_border_canvas = Image.new("RGBA", AB_circled.size,color=(255,)*4)
    inner_border_canvas = Image.new("RGBA", AB_circled.size,color=(255,)*4)
    drawCircleBorder(inner_border_canvas, inner_bbox, 1)
    drawCircleBorder(outer_border_canvas, outer_bbox, 1)

    tmp = Image.composite(AB_circled, inner_border_canvas, inner_border_canvas)
    return Image.composite(tmp, outer_border_canvas, outer_border_canvas)

def writeTagPlusValue(canvas, position, string_tag, value):
    """
    Writes a numeric value into an image, preceded with a bold tag.

    @param canvas: PIL image to write into.
    @param position: Position where it starts to write.
    @param string_tag: The tag.
    @param value: The value we want to write.
    """
    string_value = ""
    if isinstance(value,(int,long)) :
        string_value = " %d"%value
    else:
        string_value = " %.3f"%value
    writeTagPlusStringValue(canvas, position, string_tag, string_value)

def writeTagPlusStringValue(canvas, position, string_tag, string_value):
    """
    Writes text into an image, preceded with a bold tag.

    @param canvas: PIL image to write into.
    @param position: Position where it starts to write.
    @param string_tag: The tag.
    @param string_value: The string we want to write.
    """
    draw = ImageDraw.Draw(canvas)
    width, height = draw.textsize(string_tag) #@UnusedVariable
    tag_position = position
    value_position  = (position[0]+width+10,position[1])
    draw.text(tag_position, string_tag, fill = (50,50,50), font = ImageFont.truetype("/usr/share/fonts/truetype/msttcorefonts/georgiab.ttf", 11))
    draw.text(value_position, string_value, fill = (0,0,0), font = ImageFont.truetype("/usr/share/fonts/truetype/msttcorefonts/cour.ttf", 11))

def plotDataCards(image, size, number_of_cards, max_radius, h_ball_separation, v_ball_separation, data, key_exceptions):
    """
    Writes the small data cards under each cluster sphere representation.

    @param image: Canvas to draw into (PIL image).
    @param size: Tuple with the canvas size
    @param number_of_cards: Number of cards we have to write (same as number of clusters)
    @param max_radius: Is the radius of the biggest possible sphere.
    @param h_ball_separation: Separation of each sphere in the x axis.
    @param v_ball_separation: Separation of each sphere in the y axis (usually enough to fit the 'data card')
    @param data: A dictionary containing all the data with the form {data_key:(data_label, data_array)},
    for example: {"cluster_sizes":("size", [4,6,2,5,3,8])}
    @param key_exceptions: Those keys of the data dictionary we want to skip.

    @return: A PIL image with the written data cards.

    """
    w, h = size
    max_number_of_cards_in_a_row = int(w / (float(max_radius)+h_ball_separation))
    max_number_of_cards_in_a_column = int(h / (float(max_radius)+v_ball_separation))
    draw = ImageDraw.Draw(image)
    for i in range(max_number_of_cards_in_a_column):
        for j in range(max_number_of_cards_in_a_row):
            k = i*max_number_of_cards_in_a_row+j
            if k < number_of_cards:
                offset = 0
                for key in data:
                    if not key in key_exceptions:
                        datum = data[key]
                        desc = datum[0]
                        value = datum[1][k]
                        text_to_draw = desc+" %.3f"%value
                        width, height = draw.textsize(text_to_draw) #@UnusedVariable
                        offset += height + 3
                        position = (j*(max_radius+h_ball_separation), i*(max_radius+v_ball_separation)+max_radius+offset)
                        writeTagPlusValue(image, position, desc, value)
    return image

def plotSummaryTable(size, clustering_statistics_dic, per_cluster_statistics, total_elements, trajectory_comparison):
    """
    Writes the summary table, used in both the simple and extended output.

    @param size: Size of the summary table.
    @param clustering_statistics_dic:
    @param per_cluster_statistics:
    @param total_elements:
    @param trajectory_comparison: If true, it means we are performing a trajectory comparison, and the trajectory comparison
    specific data is added.

    @return: A PIL image with the written summary table.
    """
    canvas = Image.new("RGBA", size,color=(255,)*4)
    draw = ImageDraw.Draw(canvas)
    lines = []
    lines.append(("Num. elements:","   %d"%(total_elements)))
    if trajectory_comparison:

        lines.append(("Num. elements A:","   %d (%.2f%%, self: %.2f%%)"%(clustering_statistics_dic["number_elements_pure_A"],\
                                                                     clustering_statistics_dic["elems_percent_pure_A"],\
                                                                     clustering_statistics_dic["elems_self_percent_pure_A"])))
        lines.append(("Num. elements B:","   %d (%.2f%%, self: %.2f%%)"%(clustering_statistics_dic["number_elements_pure_B"],\
                                                                     clustering_statistics_dic["elems_percent_pure_B"],\
                                                                     clustering_statistics_dic["elems_self_percent_pure_B"])))
        lines.append(("Num. elements Mixed:","   %d (%.2f%%)"%(clustering_statistics_dic["number_elements_mixed"],\
                                                                     clustering_statistics_dic["elems_percent_mixed"])))

    lines.append(("Mean cluster size:","   %.2f (%.2f)"%(numpy.mean(per_cluster_statistics["cluster_sizes"][1]),\
                                                   numpy.std(per_cluster_statistics["cluster_sizes"][1]))))

    lines.append(("Mean distance to center:","   %.2f (%.2f)"%(numpy.mean(per_cluster_statistics["cluster_mean_distances"][1]),\
                                                   numpy.std(per_cluster_statistics["cluster_mean_distances"][1]))))

    lines.append(("Mean max. distance:","   %.2f (%.2f)"%(numpy.mean(per_cluster_statistics["cluster_max_distances"][1]),\
                                                   numpy.std(per_cluster_statistics["cluster_max_distances"][1]))))

    lines.append(("Mean dispersion:","   %.2f (%.2f)"%(numpy.mean(per_cluster_statistics["cluster_dispersions"][1]),\
                                                   numpy.std(per_cluster_statistics["cluster_dispersions"][1]))))
    if trajectory_comparison:
        lines.append(("Mean center differences:","   %.2f (%.2f)"%(numpy.mean(remove_zeros(per_cluster_statistics["center_differences"][1])),\
                                                       numpy.std(remove_zeros(per_cluster_statistics["center_differences"][1])))))
    initial_pos = (10,30)
    offset = 0
    for l in lines:
        position = (initial_pos[0],initial_pos[1]+offset)
        writeTagPlusStringValue(canvas, position, l[0], l[1])
        offset += draw.textsize(l[1])[1]+5
    return canvas

def fig2data (fig):
    """
    Convert a Matplotlib figure to a 4D numpy array with RGBA channels and return it

    @param fig: a matplotlib figure

    @return: a numpy 3D array of RGBA values
    """
    # draw the renderer
    fig.canvas.draw ( )

    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    buf = numpy.fromstring ( fig.canvas.tostring_argb(), dtype=numpy.uint8 )
    buf.shape = ( w, h, 4 )

    # canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
    buf = numpy.roll ( buf, 3, axis = 2 )
    return buf

def fig2img (fig):
    """
    Convert a Matplotlib figure to a PIL Image in RGBA format and return it

    @param fig: a matplotlib figure

    @return: a Python Imaging Library ( PIL ) image
    """
    # put the figure pixmap into a numpy array
    buf = fig2data ( fig )
    w, h, d = buf.shape #@UnusedVariable
    return Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )

def normalize(mylist, max_value):
    """
    Rescales the values of an array to be in the range (-max_value, max_value).

    @param mylist: The array to normalize.
    @param max_value: Maximum value in the range.

    @return: The normalized list (numpy.array).
    """
    normlist = []
    max_in_list = numpy.max(mylist)
    for d in mylist:
        normlist.append(d*float(max_value)/float(max_in_list))
    return normlist

def normalize_in_range(mylist, min_value, max_value):
    """
    Rescales all the values of an array to be in the range (min_value, max_value).

    @param mylist: The array to normalize.
    @param min_value: Minimum value we want to have in the list.
    @param max_value: Maximum value in the range.

    @return: The normalized list.
    """
    normlist = []
    max_in_list = numpy.max(mylist)
    effective_range = max_value - min_value
    for d in mylist:
        normlist.append(min_value + effective_range*(d/float(max_in_list)))
    return normlist

def remove_zeros(mylist):
    """
    Creates a new list which is equal to the input list but without zeros.

    @param mylist: The array to delete zeros.

    @return: The list without zeros.
    """
    myretlist = []
    for elem in mylist:
        if elem!=0:
            myretlist.append(elem)
    return myretlist

def tuple_to_int(t):
    """
    Creates a copy of a tuple, populated with int casts of its elements.

    @param t: The tuple to cast.

    @return: A copy of the tuple with int elements.
    """
    mytmplist = []
    for element in t:
        mytmplist.append(int(element))
    return tuple(mytmplist)

def shorten_name(name, max_length = 10):
    """
    Makes a string shorter, leaving only certain quantity of the last characters, preceded by '...'.

    @param name: The string to shorten.
    @param max_length: Maximum length of the resulting string.

    @return: A string with max_lenght characters plus '...'
    """
    if len(name) > max_length:
        return  "..."+name[-max_length:]
    else:
        return name

