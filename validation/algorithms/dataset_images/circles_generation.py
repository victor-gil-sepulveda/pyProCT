import numpy as np

# from http://wersdoerfer.com/~jochen/s9y/index.php?/archives/109-spectral-clustering-with-python.html

def get_noise(stddev=0.25, numpoints=150):  
    # 2d gaussian random noise  
    x = np.random.normal(0, stddev, numpoints)  
    y = np.random.normal(0, stddev, numpoints)  
    return np.column_stack((x, y))  
  
def get_circle(center=(0.0, 0.0), r=1.0, numpoints=150):  
    # use polar coordinates to get uniformly distributed points  
    step = np.pi * 2.0 / numpoints  
    t = np.arange(0, np.pi * 2.0, step)  
    x = center[0] + r * np.cos(t)  
    y = center[1] + r * np.sin(t)  
    return np.column_stack((x, y))  

def circle_samples():  
    circles = []  
    for radius in (1.0, 2.8, 5.0):  
        circles.append(get_circle(r=radius) + get_noise())  
    return np.vstack(circles) 
    
points = circle_samples()
for p in points:
    print "%.3f %.3f"%(p[0], p[1])
