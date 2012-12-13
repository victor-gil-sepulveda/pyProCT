'''
Created on 15/10/2012

@author: victor
'''

import numpy
import Image, ImageDraw
import matplotlib.pyplot as plt
import pylab
from pyproclust.tools.commonTools import list2ListWoZeros, normalize,\
    normalizeInRange
import ImageFont
from pyRMSD.condensedMatrix import CondensedMatrix

def matrixToImage(condensed_distance_matrix, matrix_image_file, max_size = (2048,2048)):
    # Normalize
    contents = condensed_distance_matrix.get_data()
    _max = numpy.max(contents)
    _min = numpy.min(contents)
    
    norm_contents = (contents - _min) / (_max - _min)
    
    norm_condensed = CondensedMatrix(norm_contents)
    _max = numpy.max(norm_contents)
    _min = numpy.min(norm_contents)
    
    complete = numpy.zeros([norm_condensed.row_length,norm_condensed.row_length] ,dtype=numpy.float)
    
    for i in range(norm_condensed.row_length-1):
        for j in range(i+1,norm_condensed.row_length):
            complete[i][j] = norm_condensed[i,j]
            complete[j][i] = norm_condensed[i,j]

    fig = plt.figure()
    plt.gray()
    ax = fig.add_subplot(111)
    ax.imshow(complete, interpolation='nearest')
    plt.savefig(matrix_image_file)
    
def shortenName(name):
    if len(name) > 10:
        return  "..."+name[-10:] 
    else:
        return name
    
def pieChartCreation(graph_size, fracs, name1, name2, colors):
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
    pylab.title(shortenName(name1)+" vs "+shortenName(name2))
    return fig2img(fig)

def barGraphCreation(A_sizes, B_sizes, cluster_sizes, types, total_size, colors, graph_size):
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

def tuple2Int(t):
    mytmplist = []
    for element in t:
        mytmplist.append(int(element))
    return tuple(mytmplist)
    
def ballGraph(number_of_balls, size, max_box_diameter, h_ball_separation, v_ball_separation, sizes, alphas, colors, cluster_types, A_sizes, B_sizes):
    w,h = size
    max_number_of_balls_in_a_row  = 0
    max_number_of_balls_in_a_column = 0
    max_tmp_diameter = max_box_diameter - 4
    while  max_number_of_balls_in_a_row * max_number_of_balls_in_a_column <= number_of_balls:
        max_number_of_balls_in_a_row = int(w / (float(max_tmp_diameter)+h_ball_separation))
        max_number_of_balls_in_a_column = int(h / (float(max_tmp_diameter)+v_ball_separation))
        max_tmp_diameter -= 1
    
    # Normalize ball sizes and alphas
    norm_sizes = normalize(sizes,max_box_diameter)
    # Reverse values and normalize between 1 and 0.3 (0.3 is the alpha of the cluster with most dispersion)
    norm_alphas = normalizeInRange(1-numpy.array(normalizeInRange(alphas,0,1)),0.3,1)
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
                    canvas.paste(ball,tuple2Int(pos))
                else:
                    ball_with_inner_circles = circledPercentages(ball, (centering_offset,centering_offset), this_diameter, 0.3, A_sizes[k], B_sizes[k], colors)
                    canvas.paste(ball_with_inner_circles,tuple2Int(pos))
    return canvas,max_tmp_diameter

def drawBall(canvas_size, ball_diameter, fill_color, alpha):
    canvas = Image.new("RGBA", canvas_size,color=(255,)*4)
    draw = ImageDraw.Draw(canvas)
    centering_offset = int(max(canvas_size)/2. - ball_diameter/2)
    position = (centering_offset,centering_offset,centering_offset+ball_diameter,centering_offset+ball_diameter)
    draw.ellipse(position,fill=fill_color)
    blending_canvas = Image.new("RGBA", canvas_size,color=(255,)*4)
    border_canvas = Image.new("RGBA", canvas_size,color=(255,)*4)
    drawCircleBorder(border_canvas, position, 2)
    return Image.composite(Image.blend(blending_canvas,canvas,alpha), border_canvas, border_canvas)

def drawCircleBorder(canvas,bbox,width):
    draw = ImageDraw.Draw(canvas)
    border_position = (bbox[0]-width,bbox[1]-width,bbox[2]+width,bbox[3]+width)
    draw.ellipse(border_position,fill=(0,)*4)
    draw.ellipse(tuple2Int(bbox),fill=(255,)*4)
    
def circledPercentages(ball_canvas, position, ball_diameter, radius_percent, sizeA, sizeB, colors):
    
    ball_radius = ball_diameter/2.
    width = radius_percent * ball_radius
    
    inner_bbox = (position[0]+0.5*ball_radius,position[1]+0.5*ball_radius,\
                  position[0]+ball_diameter-(0.5*ball_radius),position[1]+ball_diameter-(0.5*ball_radius))
    
    outer_bbox = (inner_bbox[0]-width, inner_bbox[1]-width,\
                  inner_bbox[2]+width, inner_bbox[3]+width)
    
    A_mask = Image.new("RGBA", ball_canvas.size,color=(255,)*4)
    drawCircleBorder(A_mask,tuple2Int(inner_bbox),width)
    A_color = Image.new("RGBA", ball_canvas.size,color=colors['A'])
    A_circled =  Image.composite(ball_canvas, A_color, A_mask)
    
    B_mask = Image.new("RGBA", ball_canvas.size,color=(255,)*4)
    drawCircleBorder(B_mask,tuple2Int(inner_bbox),width)
    B_color = Image.new("RGBA", ball_canvas.size,color=colors['B'])
    draw = ImageDraw.Draw(B_mask)
    B_percent = sizeB / float(sizeA+sizeB)
    draw.pieslice(tuple2Int(outer_bbox), 0, int(360 *(1- B_percent)), fill = (255,)*4)
    
    # Separation line from A-B
    draw = ImageDraw.Draw(B_color)
    draw.pieslice(tuple2Int(outer_bbox), -1, int(360 *(1- B_percent))+1, outline = (0,)*4 )
    
    AB_circled = Image.composite(A_circled, B_color, B_mask)
    outer_border_canvas = Image.new("RGBA", AB_circled.size,color=(255,)*4)
    inner_border_canvas = Image.new("RGBA", AB_circled.size,color=(255,)*4)
    drawCircleBorder(inner_border_canvas, inner_bbox, 1)
    drawCircleBorder(outer_border_canvas, outer_bbox, 1)
    
    tmp = Image.composite(AB_circled, inner_border_canvas, inner_border_canvas)
    return Image.composite(tmp, outer_border_canvas, outer_border_canvas)

def writeTagPlusValue(canvas, position, string_tag, value):
    string_value = ""
    if isinstance(value,(int,long)) :
        string_value = " %d"%value
    else:
        string_value = " %.3f"%value
    writeTagPlusStringValue(canvas, position, string_tag, string_value)

def writeTagPlusStringValue(canvas, position, string_tag, string_value):
    draw = ImageDraw.Draw(canvas)
    width, height = draw.textsize(string_tag) #@UnusedVariable
    tag_position = position
    value_position  = (position[0]+width+10,position[1])
    draw.text(tag_position, string_tag, fill = (50,50,50), font = ImageFont.truetype("/usr/share/fonts/truetype/msttcorefonts/georgiab.ttf", 11))
    draw.text(value_position, string_value, fill = (0,0,0), font = ImageFont.truetype("/usr/share/fonts/truetype/msttcorefonts/cour.ttf", 11))
    
def plotDataCards(image, size, number_of_cards, max_radius, h_ball_separation, v_ball_separation, data, key_exceptions):
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
        lines.append(("Mean center differences:","   %.2f (%.2f)"%(numpy.mean(list2ListWoZeros(per_cluster_statistics["center_differences"][1])),\
                                                       numpy.std(list2ListWoZeros(per_cluster_statistics["center_differences"][1])))))
    initial_pos = (10,30)
    offset = 0
    for l in lines:
        position = (initial_pos[0],initial_pos[1]+offset)
        writeTagPlusStringValue(canvas, position, l[0], l[1])
        offset += draw.textsize(l[1])[1]+5
    return canvas

def fig2data ( fig ):
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
    
def fig2img ( fig ):
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
    normlist = []
    max_in_list = numpy.max(mylist)
    for d in mylist:
        normlist.append(d*float(max_value)/float(max_in_list))
    return normlist

def normalizeInRange(mylist, min_value, max_value):
    normlist = []
    max_in_list = numpy.max(mylist)
    effective_range = max_value - min_value
    for d in mylist:
        normlist.append(min_value + effective_range*(d/float(max_in_list)))
    return normlist

def list2ListWoZeros(mylist):
    myretlist = []
    for elem in mylist:
        if elem!=0:
            myretlist.append(elem)
    return myretlist

if __name__ == '__main__':
    image, max_radius = ballGraph(150, (400,400), 150, 10, 40, range(150),range(150))
    image.show()
    data = {"lol":("lol:",range(150)),"lel":("lel:",range(150,300))}
    
    plotDataCards(image, (400,400), 150, max_radius, 10, 40, data)
    image.show()
    