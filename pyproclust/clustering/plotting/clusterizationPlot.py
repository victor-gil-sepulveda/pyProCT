'''
Created on 23/03/2012

@author: victor
'''

import scipy.cluster.hierarchy as hcluster
import matplotlib.pylab

def plot_dendogram(self,link_matrix,save_path = None):
    """
    Plots a dendogram using matplotlib and required writes it to save_path.
    Very difficult to test :/
    """
    fig = matplotlib.pylab.figure(figsize=(10, 8), dpi=160)
    hcluster.dendrogram(link_matrix)
    matplotlib.pylab.show()
    if save_path:
        fig.savefig(save_path,format="png")      

