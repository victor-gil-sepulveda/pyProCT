'''
Created on 1/12/2014

@author: victor
'''
import math

def sq_distance(a,b):
    return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2

def angular_rmsd(a,b):
    return math.sqrt((angular_increment(a-b)**2).sum()/len(a))

def to_0_2PI_range(angle):
    """
    Converts an angle from the  [-pi, pi] range to the [0, 2pi] range.
    :param angle: The angle to be converted.
    
    :return: The converted angle in [-pi, pi] range.
    """
    if angle >= 0 :
        return angle
    else:
        return (2*math.pi)+angle;

def angular_increment(angle1, angle2):
    return min(angle1-angle2, to_0_2PI_range(angle1)-to_0_2PI_range(angle2))

def angular_distance(angle1, angle2):
    """
    angle1, angle2 , angles in -pi, pi range. 
    
    Returns the angular increment defined by these angles.
    """
    return min(abs(angle1-angle2), abs(to_0_2PI_range(angle1)-to_0_2PI_range(angle2)))

def calc_dihedral(a1_coords, a2_coords, a3_coords, a4_coords):
    xa, ya, za = a1_coords
    xb, yb, zb = a2_coords
    xc, yc, zc = a3_coords
    xd, yd, zd = a4_coords
    
    v1 = xc-xd
    v2 = yb-yc
    v3 = xb-xc
    v4 = yc-yd
    v5 = v3*v4-v1*v2
    v6 = ya-yb
    v7 = xa-xb
    v8 = v7*v2-v3*v6
    v9 = za-zb
    v10 = zb-zc
    v11 = v3*v9-v7*v10
    v12 = v6*v10-v2*v9
    v13 = 1/math.sqrt(v8*v8+v12*v12+v11*v11)
    v14 = 1/math.sqrt(v3*v3+v2*v2+v10*v10)
    v15 = zc-zd;
    v16 = v1*v10-v3*v15
    v17 = v2*v15-v4*v10
    v18 = 1/math.sqrt(pow(v5,2)+pow(v17,2)+pow(v16,2))

    return math.atan2((v13*v14*v11*v10-v8*v2*v13*v14)*v18*v17+(v3*v8*v13*v14-v13*v14*v12*v10)*v18*v16+v5*(v2*v13*v14*v12-v3*v13*v14*v11)*v18,v13*v12*v18*v17+v13*v11*v18*v16+v8*v5*v13*v18)
     