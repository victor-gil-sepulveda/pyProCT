import numpy as np

# from http://wersdoerfer.com/~jochen/s9y/index.php?/archives/109-spectral-clustering-with-python.html

# def get_noise(stddev=0.25, numpoints=150):  
#     # 2d gaussian random noise  
#     x = np.random.normal(0, stddev, numpoints)  
#     y = np.random.normal(0, stddev, numpoints)  
#     return np.column_stack((x, y))  
#   
# def get_circle(center=(0.0, 0.0), r=1.0, numpoints=150):  
#     # use polar coordinates to get uniformly distributed points  
#     step = np.pi * 2.0 / numpoints  
#     t = np.arange(0, np.pi * 2.0, step)  
#     x = center[0] + r * np.cos(t)  
#     y = center[1] + r * np.sin(t)  
#     return np.column_stack((x, y))  
# 
# def circle_samples():  
#     circles = []  
#     for radius in (1.0, 2.8, 5.0):  
#         circles.append(get_circle(r=radius) + get_noise())  
#     return np.vstack(circles) 
#     
# points = circle_samples()
# for p in points:
#     print "%.3f %.3f"%(p[0], p[1])

import sys

import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import eig
from scipy.cluster.vq import kmeans2
from scipy.sparse.linalg import eigen
from scipy.spatial.kdtree import KDTree

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

def radial_kernel(c=1.5):
    def inner(a, b):
        d = a - b
        return np.exp((-1 * (np.sqrt(np.dot(d, d.conj()))**2)) / c)
    return inner

def circle_samples():
    circles = []
    for radius in (1.0, 2.8, 5.0):
        circles.append(get_circle(r=radius) + get_noise())
    return np.vstack(circles)

def mutual_knn(points, n=10, distance=radial_kernel()):
    knn = {}
    kt = KDTree(points)
    for i, point in enumerate(points):
        # cannot use euclidean distance directly
        for neighbour in kt.query(point, n + 1)[1]:
            if i != neighbour:
                knn.setdefault(i, []).append(
                    (distance(point, points[neighbour]), neighbour))
    return knn

def get_distance_matrix(knn):
    n = len(knn)
    W = np.zeros((n, n))
    for point, nearest_neighbours in knn.iteritems():
        for distance, neighbour in nearest_neighbours:
            W[point][neighbour] = distance
    return W

def rename_clusters(idx):
    # so that first cluster has index 0
    num = -1
    seen = {}
    newidx = []
    for id in idx:
        if id not in seen:
            num += 1
            seen[id] = num
        newidx.append(seen[id])
    return np.array(newidx)

def cluster_points(L):
    # sparse eigen is a little bit faster than eig
    #evals, evcts = eigen(L, k=15, which="SM")
    evals, evcts = eig(L)
    evals, evcts = evals.real, evcts.real
    edict = dict(zip(evals, evcts.transpose()))
    evals = sorted(edict.keys())
    # second and third smallest eigenvalue + vector
    Y = np.array([edict[k] for k in evals[1:6]]).transpose()
    res, idx = kmeans2(Y, 6, minit='random')
    return evals[:15], Y, rename_clusters(idx)

def change_tick_fontsize(ax, size):
    for tl in ax.get_xticklabels():
        tl.set_fontsize(size)
    for tl in ax.get_yticklabels():
        tl.set_fontsize(size)

def get_colormap():
    # map cluster label to color (0, 1, 2) -> (orange, blue, green)
    from matplotlib.colors import ListedColormap
    orange = (0.918, 0.545, 0.0)
    blue = (0.169, 0.651, 0.914)
    green = (0.0, 0.58, 0.365)
    return ListedColormap([orange, blue, green])

def plot_circles(ax, points, idx, colormap):
    plt.scatter(points[:,0], points[:,1], s=10, c=idx, cmap=colormap,
        alpha=0.9, facecolors="none")
    plt.xlabel("x1", fontsize=8)
    plt.ylabel("x2", fontsize=8)
    change_tick_fontsize(ax, 8 )
    plt.ylim(-6, 6)
    plt.xlim(-6, 6)

def plot_eigenvalues(ax, evals):
    plt.scatter(np.arange(0, len(evals)), evals,
        c=(0.0, 0.58, 0.365), linewidth=0)
    plt.xlabel("Number", fontsize=8)
    plt.ylabel("Eigenvalue", fontsize=8)
    plt.axhline(0, ls="--", c="k")
    change_tick_fontsize(ax, 8 )

def plot_eigenvectors(ax, Y, idx, colormap):
    from matplotlib.ticker import MaxNLocator
    from mpl_toolkits.axes_grid import make_axes_locatable
    divider = make_axes_locatable(ax)
    ax2 = divider.new_vertical(size="100%", pad=0.05)
    fig1 = ax.get_figure()
    fig1.add_axes(ax2)
    ax2.set_title("Eigenvectors", fontsize=10)
    ax2.scatter(np.arange(0, len(Y)), Y[:,0], s=10, c=idx, cmap=colormap,
        alpha=0.9, facecolors="none")
    ax2.axhline(0, ls="--", c="k")
    ax2.yaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.axhline(0, ls="--", c="k")
    ax.scatter(np.arange(0, len(Y)), Y[:,1], s=10, c=idx, cmap=colormap,
        alpha=0.9, facecolors="none")
    ax.set_xlabel("index", fontsize=8)
    ax2.set_ylabel("2nd Smallest", fontsize=8)
    ax.set_ylabel("3nd Smallest", fontsize=8)
    change_tick_fontsize(ax, 8 )
    change_tick_fontsize(ax2, 8 )
    for tl in ax2.get_xticklabels():
        tl.set_visible(False)

def plot_spec_clustering(ax, Y, idx, colormap):
    plt.title("Spectral Clustering", fontsize=10)
    plt.scatter(Y[:,0], Y[:,1], c=idx, cmap=colormap, s=10, alpha=0.9,
        facecolors="none")
    plt.xlabel("Second Smallest Eigenvector", fontsize=8)
    plt.ylabel("Third Smallest Eigenvector", fontsize=8)
    change_tick_fontsize(ax, 8 )

def plot_figure(points, evals, Y, idx):
    colormap = get_colormap()
    fig = plt.figure(figsize=(6, 5.5))

    fig.subplots_adjust(wspace=0.4, hspace=0.3)
    ax = fig.add_subplot(2, 2, 1)
    plot_circles(ax, points, idx, colormap)

    ax = fig.add_subplot(2, 2, 2)
    plot_eigenvalues(ax, evals)

    ax = fig.add_subplot(2, 2, 3)
    plot_eigenvectors(ax, Y, idx, colormap)

    ax = fig.add_subplot(2, 2, 4)
    plot_spec_clustering(ax, Y, idx, colormap)

    plt.show()

def main(args):
    points = circle_samples()
    knn_points = mutual_knn(points)
    W = get_distance_matrix(knn_points)
    G = np.diag([sum(Wi) for Wi in W])

    # unnormalized graph Laplacian
    L = G - W
    evals, Y, idx = cluster_points(L)

    plot_figure(points, evals, Y, idx)

if __name__ == "__main__":
    main(sys.argv)
