import sys
import json
from pyproct.tools.commonTools import convert_to_utf8
from pyproct.clustering.clustering import Clustering
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.spines import Spine
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection

N = 9

################################################################
# From http://matplotlib.org/examples/api/radar_chart.html
################################################################
def unit_poly_verts(theta):
    """Return vertices of polygon for subplot axes.

    This polygon is circumscribed by a unit circle centered at (0.5, 0.5)
    """
    x0, y0, r = [0.5] * 3
    verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
    return verts

def radar_factory(num_vars, frame='circle'):
    """Create a radar chart with `num_vars` axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle' | 'polygon'}
        Shape of frame surrounding axes.

    """
    # calculate evenly-spaced axis angles
    theta = 2*np.pi * np.linspace(0, 1-1./num_vars, num_vars)
    # rotate theta such that the first axis is at the top
    theta += np.pi/2

    def draw_poly_patch(self):
        verts = unit_poly_verts(theta)
        return plt.Polygon(verts, closed=True, edgecolor='k')

    def draw_circle_patch(self):
        # unit circle centered on (0.5, 0.5)
        return plt.Circle((0.5, 0.5), 0.5)

    patch_dict = {'polygon': draw_poly_patch, 'circle': draw_circle_patch}
    if frame not in patch_dict:
        raise ValueError('unknown value for `frame`: %s' % frame)

    class RadarAxes(PolarAxes):

        name = 'radar'
        # use 1 line segment to connect specified points
        RESOLUTION = 1
        # define draw_frame method
        draw_patch = patch_dict[frame]

        def fill(self, *args, **kwargs):
            """Override fill so that line is closed by default"""
            closed = kwargs.pop('closed', True)
            return super(RadarAxes, self).fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super(RadarAxes, self).plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.concatenate((x, [x[0]]))
                y = np.concatenate((y, [y[0]]))
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(theta * 180/np.pi, labels)

        def _gen_axes_patch(self):
            return self.draw_patch()

        def _gen_axes_spines(self):
            if frame == 'circle':
                return PolarAxes._gen_axes_spines(self)
            # The following is a hack to get the spines (i.e. the axes frame)
            # to draw correctly for a polygon frame.

            # spine_type must be 'left', 'right', 'top', 'bottom', or `circle`.
            spine_type = 'circle'
            verts = unit_poly_verts(theta)
            # close off polygon by repeating first vertex
            verts.append(verts[0])
            path = Path(verts)

            spine = Spine(self, spine_type, path)
            spine.set_transform(self.transAxes)
            return {'polar': spine}
    register_projection(RadarAxes)
    return theta
################################################################


def gen_data(cluster, index_to_interpolation):
    cluster.percents = {}
    for element in cluster.all_elements:
        try:
            cluster.percents[index_to_interpolation[element]] += 1
            a,b = index_to_interpolation[element]
            cluster.percents[(b,a)] += 1
        except KeyError:
            cluster.percents[index_to_interpolation[element]] = 1
            a,b = index_to_interpolation[element]
            cluster.percents[(b,a)] = 1

    for percent in cluster.percents:
        cluster.percents[percent] /= 20.

    data = {}
    data["column names"] = [str(x) for x in range(N)]
    for i in range(0,N):
        data[str(i)] = []
        for j in range(0,N):
            if (i,j) in cluster.percents:
                data[str(i)].append( cluster.percents[(i,j)])
            else:
                data[str(i)].append(0)

    return data



if __name__ == '__main__':
    results = convert_to_utf8(json.loads(open(sys.argv[1]).read()))
    best_clustering_id =results["best_clustering"]
    best_clustering_dic = results["selected"][best_clustering_id]
    num_clusters = best_clustering_dic["clustering"]["number_of_clusters"]
    clustering = Clustering.from_dic(best_clustering_dic["clustering"])
    file_frames = int(sys.argv[2])

    # generate a map element -> interpolation
    index_to_interpolation = {}
    acc = 0
    for i in range(0, file_frames-1):
        for j in range(i+1, file_frames):
            for k in range(20):
                index_to_interpolation[acc] = (i,j)
                acc += 1


    for cluster in clustering.clusters:
        colors = iter(cm.rainbow(np.linspace(0, 1, N)))
        theta = radar_factory(N, frame='polygon')

        data = gen_data(cluster, index_to_interpolation)
        spoke_labels = data.pop('column names')

        fig = plt.figure(figsize=(9, 9))
        fig.subplots_adjust(wspace = 0.25, hspace = 0.20, top = 0.85, bottom = 0.05)

        # Plot the four cases from the example data on separate axes
        for n, title in enumerate(data.keys()):
            ax = fig.add_subplot(3, 3, n+1, projection='radar')
            ax.set_rmax(1.)
            ax.autoscale(False)
            plt.rgrids([0.5])
            ax.set_title(title, weight='bold', size='medium', position=(0.5, 1.1),
                         horizontalalignment='center', verticalalignment='center')
            color = next(colors)
            sorted_data_indices =  np.argsort(data[title])
            sorted_data = np.array(data[title])[sorted_data_indices]
            sorted_labels = np.array(spoke_labels)[sorted_data_indices]
            ax.plot(theta, sorted_data, color='k')
            ax.fill(theta, sorted_data, facecolor=color, alpha=0.25)
            ax.set_varlabels(sorted_labels)

        plt.figtext(0.5, 0.965, cluster.id, ha='center', color='black', weight='bold', size='large')
        plt.show()
