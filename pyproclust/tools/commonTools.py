'''
Created on 16/03/2012

@author: victor
'''
import sys
import numpy

def merge_files(file_handler_list, merged_handler, verbose = True):
    """
    Merges the files which are in the file_handler_list to write them, line per line,
    into the merged_handler.
    If the files come from the generate_rmsd_matrix_in_parallel_with_vmd function, 
    they have to be sorted in order to recreate a correct matrix.
    """
    total_files = len(file_handler_list)
    current_file = 1
    if verbose:
        print ""
    for f in file_handler_list:
        if verbose:
            print "Processing file",current_file,"of",total_files
        for line in f:
            merged_handler.write(line)
        current_file = current_file +1 

def open_file_handler_list(files):
    """
    Opens a list of files and returns its handlers.
    """
    handlers = []
    for f in files:
        handlers.append(open(f,"r"))
    return handlers

def print_and_flush(this_string):
    """
    No comments...
    """
    sys.stdout.write(this_string)
    sys.stdout.flush()
    
def vararg_callback(option, opt_str, value, parser):
    """
    To use when parsing command line parameters (parses a list of float numbers)
    """
    assert value is None
    value = []
    
    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False
    
    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(float(arg))
    
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)

def gen_consecutive_ranges(num_elems_1,num_elems_2):
    """
    Generates two consecutive ranges from 0 to num_elems_1-1 and from num_elems_1 to
    num_elems_1+num_elems_2 -1 . 
    """
    return range(num_elems_1),range(num_elems_1,num_elems_1+num_elems_2)

def list2ListWoZeros(mylist):
    myretlist = []
    for elem in mylist:
        if elem!=0:
            myretlist.append(elem)
    return myretlist

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

# from http://mjijackson.com/2008/02/rgb-to-hsl-and-rgb-to-hsv-color-model-conversion-algorithms-in-javascript
#/**
# * Converts an RGB color value to HSL. Conversion formula
# * adapted from http://en.wikipedia.org/wiki/HSL_color_space.
# * Assumes r, g, and b are contained in the set [0, 255] and
# * returns h, s, and l in the set [0, 1].
# *
# * @param   Number  r       The red color value
# * @param   Number  g       The green color value
# * @param   Number  b       The blue color value
# * @return  Array           The HSL representation
# */
def rgbToHsl(r_unorm, g_unorm, b_unorm):
    r = r_unorm/255
    g = g_unorm/255
    b = b_unorm/255
    maxc = max([r, g, b])
    minc = min([r, g, b])
    
    h = (maxc + minc) / 2.
    s = (maxc + minc) / 2.
    l = (maxc + minc) / 2.
    
    if maxc == minc:
        h = 0
        s = 0 # achromatic
    else:
        d = maxc - minc
        if l > 0.5:
            s = d / (2 - maxc - minc)
        else:
            s = d / (maxc + minc)
        if maxc == r: 
            modif = 0
            if g < b:
                modif = 6
            h = (g - b) / d + modif
        elif maxc == g: 
            h = (b - r) / d + 2
        elif maxc == b:
            h = (r - g) / d + 4
        h /= 6;
    return (h, s, l);


#/**
# * Converts an HSL color value to RGB. Conversion formula
# * adapted from http://en.wikipedia.org/wiki/HSL_color_space.
# * Assumes h, s, and l are contained in the set [0, 1] and
# * returns r, g, and b in the set [0, 255].
# *
# * @param   Number  h       The hue
# * @param   Number  s       The saturation
# * @param   Number  l       The lightness
# * @return  Array           The RGB representation
# */
def hslToRgb(h, s, l):
    r = 0
    g = 0
    b = 0
    
    if s == 0 :
        r = g = b = l  # achromatic
    else:
        def hue2rgb(p, q, t):
            if t < 0: t += 1
            if t > 1: t -= 1
            if t < 1./6: return p + (q - p) * 6 * t
            if t < 1./2: return q
            if t < 2./3: return p + (q - p) * (2/3 - t) * 6
            return p
        
        if l < 0.5:
            q = l * (1 + s)
        else:
            q = l + s - l * s
        p = 2 * l - q
        r = hue2rgb(p, q, h + 1./3)
        g = hue2rgb(p, q, h)
        b = hue2rgb(p, q, h - 1/3)
    
    return (r * 255, g * 255, b * 255)


