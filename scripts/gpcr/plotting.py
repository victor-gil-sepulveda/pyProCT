import matplotlib
import matplotlib.pyplot as plt
import numpy
import matplotlib.cm as cm

def plot_metrics(filename, metrics, x_label, y_label):
    markersize = 10.
    fig, ax = plt.subplots()
    ax.scatter(metrics.T[0],metrics.T[1], s = markersize)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig(filename)

def calculate_rectangle(dl, ur, extra = 0, xy_ratio = 1.):
    dlx, dly = dl[0] - extra*xy_ratio, dl[1] - extra
    urx, ury = ur[0] + extra*xy_ratio, ur[1] + extra
    return {
            'center': (dlx + (urx - dlx)/2., dly + (ury - dly)/2.),
            'width': urx - dlx,
            'height':ury - dly,
            'down left corner': (dlx, dly)
            }

def plot_clusters(filename, metrics, scores, best, plot_rectangles = False):
    markersize = 10.
    fig, ax = plt.subplots()
    pixel_markersize = markersize/fig.dpi

    # get scale ratio of measures
    rmsd_range = numpy.max(metrics.T[0]) - numpy.min(metrics.T[0])
    be_range = numpy.max(metrics.T[1]) - numpy.min(metrics.T[1])
    ratio = rmsd_range/float(be_range)
    cmap = matplotlib.colors.ListedColormap ( numpy.random.rand ( 256,3))

    for i, (score, cluster) in enumerate(scores):
        filtered_metrics = metrics[cluster.all_elements]
        rectangle = calculate_rectangle((numpy.min(filtered_metrics.T[0]), numpy.min(filtered_metrics.T[1])),
                                        (numpy.max(filtered_metrics.T[0]), numpy.max(filtered_metrics.T[1])),
                                        extra = pixel_markersize,
                                        xy_ratio = ratio )
        if plot_rectangles:
            ax.add_patch(matplotlib.patches.Rectangle(rectangle["down left corner"],
                                         rectangle['width'],
                                         rectangle['height'],
                                         facecolor=cm.prism(255-i),
                                         alpha =  0.5,
                                         edgecolor = 'black',
                                         zorder=-1))

        ax.scatter(filtered_metrics.T[0],filtered_metrics.T[1], c = cmap(i) , s = markersize)

        if i == best:
            ax.annotate(r'Best %d'%i, rectangle['center'], size = 'x-small', ha='center')
        else:
            ax.annotate("%d"%i, rectangle['center'], alpha = 0.7, size = 'xx-small', ha='center')
    plt.xlabel("RMSD (A)")
    plt.ylabel("Binding Energy")
    plt.savefig(filename)
