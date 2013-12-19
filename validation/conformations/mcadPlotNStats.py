'''
Created on 19/12/2013

@author: victor
'''
import sys
from pyproct.tools.commonTools import convert_to_utf8
import json
import matplotlib.pyplot as plt
import numpy

displacement_files = {
                      "Defined Clusters":'/home/victor/Escritorio/pyProCT-Tests/Synthetic/9_natural_clusters/results/CA_displacements.json',
                      "Diffuse Clusters (1)":'/home/victor/Escritorio/pyProCT-Tests/Synthetic/9_diffuse_clusters/7/results/CA_displacements.json',
                      "Diffuse Clusters (2)":'/home/victor/Escritorio/pyProCT-Tests/Synthetic/9_diffuse_clusters/4/results/CA_displacements.json',
                      "Diffuse Clusters (3)":'/home/victor/Escritorio/pyProCT-Tests/Synthetic/9_diffuse_clusters/5/results/CA_displacements.json'
                      }
ordered_keys =["Defined Clusters", "Diffuse Clusters (1)", "Diffuse Clusters (2)", "Diffuse Clusters (3)"]

def calculate_max_min(displacements):
    all_displacements = []
    for cluster_id in displacements:
        if cluster_id != "global":
            all_displacements.append(displacements[cluster_id])
    np_displacements = numpy.array(all_displacements)
    return numpy.max(np_displacements,0),numpy.min(np_displacements,0)

if __name__ == '__main__':
    fig = plt.figure(figsize=(9, 9))
    fig.subplots_adjust(wspace = 0.25, hspace = 0.20, top = 0.85, bottom = 0.05)

    for i, file_id in enumerate(ordered_keys):
        print i, file_id
        displacements = convert_to_utf8(json.loads(open(displacement_files[file_id]).read()))
        ax = fig.add_subplot(len(displacement_files), 1, i)
        ax.set_ylim(-1,max(displacements["global"])+5)
        ax.plot(displacements["global"], color='k', linewidth =2)
        max_displacements, min_displacements = calculate_max_min(displacements)
        ax.plot(max_displacements, color='r', linewidth = 1)
        ax.plot(min_displacements, color='b', linewidth = 1)
        ax.set_title(file_id, weight='bold', size='medium', #position=(0.5, 1.1),
                         horizontalalignment='center', verticalalignment='center')
    plt.show()
