'''
Created on 17/09/2012

@author: victor
'''
from pyproclust.tools.scriptTools import get_directories_list
import re
from pylab import *
import matplotlib.pyplot as plt
import numpy

if __name__ == '__main__':
    directories = get_directories_list(".")
    kl_all = {}
    overlap_all = {}
    cluster_populations_all = {}
    element_percents = {}
    all_files = []
    all_clustering_sizes = {}
    filtered_directories = []
    
    for d in directories:
        if '_vs_' in d :
            filtered_directories.append(d)
            
    for d in filtered_directories:
        results_dir = d+"/results"
        KL_file_p =  d+"/matrix/rmsd_distrib.txt"#results_dir+"/rmsd_distrib.txt"
        overlap_file_p = results_dir+"/overlap.txt"
        try:
            KL_file = open(KL_file_p,"r")
            overlap_file = open(overlap_file_p,"r")
            
            tags = d.split("_vs_")
            traj1_name = tags[0]
            traj2_name = tags[1]
            all_files.append(traj1_name)
            all_files.append(traj2_name)
            def parse_KL_file(file_handler,name1,name2,kl_all):
                for l in file_handler:
                    if "t1 on t2" in l:
                        tags = l.split(":")
                        kl_all[name1,name2] = float(tags[1])
                    if "t2 on t1" in l:
                        tags = l.split(":")
                        kl_all[name2,name1] = float(tags[1])
                        
            parse_KL_file(KL_file,traj1_name,traj2_name,kl_all)
            KL_file.close()
            
            def parse_overlap_file(file_handler,name1,name2,overlap_all):
                
                print name1, name2
                lines = file_handler.readlines()
                
                m = re.search(r'The clustering\ has\ (\d+)\ elements\ \((\d+)\ from',lines[2])
                if m:    
                    T_elems = int(m.group(1))
                    A_elems = int(m.group(2))
                m = re.search(r'and (\d+) from',lines[2])
                if m:    
                    B_elems = int(m.group(1))
                
                A_percent = float(lines[3].split()[1][0:-1])
                B_percent = float(lines[5].split()[1][0:-1])
                mixed_percent = float(lines[7].split()[1][0:-1])

                element_percents[name1,name2] = (T_elems, A_elems, B_elems, A_percent, B_percent, mixed_percent)
                element_percents[name2,name1] = (T_elems, B_elems, A_elems, B_percent, A_percent, mixed_percent)
                
                m = re.search(r'corresponding\ to\ (\d+)\ clusters',lines[3])
                if m:    
                    A_num = int(m.group(1))
                m = re.search(r'corresponding\ to\ (\d+)\ clusters',lines[5])
                if m:    
                    B_num = int(m.group(1))
                m = re.search(r'corresponding\ to\ (\d+)\ clusters',lines[7])
                if m:    
                    mixed = int(m.group(1))
                
                
                if not name1 in all_clustering_sizes:
                    all_clustering_sizes[name1] = []
                all_clustering_sizes[name1].append(A_num+mixed)
                if not name2 in all_clustering_sizes:
                    all_clustering_sizes[name2] = []
                all_clustering_sizes[name2].append(B_num+mixed)
                
                tags = lines[9].split(":")
                mean_trans_A_over_B = float(tags[1])
                tags = lines[10].split(":")
                mean_scale_A_over_B = float(tags[1])
                
                tags = lines[12].split(":")
                mean_trans_B_over_A = float(tags[1])
                tags = lines[13].split(":")
                mean_scale_B_over_A = float(tags[1])
                    
                overlap_all[name1,name2] = (A_num,B_num,mixed,mean_trans_A_over_B,mean_scale_A_over_B)
                overlap_all[name2,name1] = (B_num,A_num,mixed,mean_trans_B_over_A,mean_scale_B_over_A)
                
                try:
                    # parse population 
                    population_line = lines[14].split(":")
                    cluster_populations = population_line[1].split()
                    
                    cluster_populations_all[name1,name2] = []
                    cluster_populations_all[name2,name1] = []
                    complementary_tags = {"A":"B","B":"A","M":"M"}
                    for cpop in cluster_populations:
                        tags = cpop.split("/")
                        if len(tags) == 3: #'M'
                            cluster_populations_all[name1,name2].append((float(tags[0])*100,float(tags[1])*100,tags[2]))
                            cluster_populations_all[name2,name1].append((float(tags[1])*100,float(tags[0])*100,complementary_tags[tags[2]]))
                        else:
                            cluster_populations_all[name1,name2].append((float(tags[0])*100,tags[1]))
                            cluster_populations_all[name2,name1].append((float(tags[0])*100,complementary_tags[tags[1]]))
                except IndexError:   
                    print "Population line not found in ",overlap_file_p
            
            parse_overlap_file(overlap_file,traj1_name,traj2_name,overlap_all)
            overlap_file.close()
            
        except IOError:
            print "file not found in ",KL_file_p,"or",overlap_file_p
    
    matplotlib.rcParams['lines.linewidth'] = 100
    # Create graphs
    def radial_graph_creation(fracs, name1, name2, trans, scale):
        all_labels = ["A","B","Mixed"]
        all_colors = ['#468C54','#4786A1','#949E5D']
        this_labels = []
        this_fracs = []
        this_colors = []
        fig = figure(figsize=(4,4))
        fig.set_facecolor('white')
        axes([0.1, 0.1, 0.8, 0.8])
        for i in range(len(fracs)):
            if fracs[i]!=0:
                this_labels.append(all_labels[i])
                this_fracs.append(fracs[i])
                this_colors.append(all_colors[i])
        pie(this_fracs, labels=this_labels, autopct='%1.1f%%', shadow=False,colors=this_colors)
        title(name1+" vs "+name2)
        
        (T_elems, A_elems, B_elems, A_percent, B_percent, mixed_percent) = element_percents[name1,name2]
        fig.text(0,0.09,"Elems: %d (A: %d B: %d)"%(T_elems,A_elems, B_elems),size = 10)
                                                                                          
        fig.text(0,0.05,"Clusters A: %d (%.2f%%), B: %d (%.2f%%), M: %d (%.2f%%) "%(fracs[0],A_percent,\
                                                                           fracs[1],B_percent,\
                                                                           fracs[2],mixed_percent),size = 10)
        fig.text(0,0.01,"KL: %0.4f Trans.: %0.3f Scale: %0.3f"%(kl_all[name1,name2],trans,scale),size = 10 )
        return fig
    
    import numpy
    def fig2data ( fig ):
        """
        @brief Convert a Matplotlib figure to a 4D numpy array with RGBA channels and return it
        @param fig a matplotlib figure
        @return a numpy 3D array of RGBA values
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
    
    import Image
     
    def fig2img ( fig ):
        """
        @brief Convert a Matplotlib figure to a PIL Image in RGBA format and return it
        @param fig a matplotlib figure
        @return a Python Imaging Library ( PIL ) image
        """
        # put the figure pixmap into a numpy array
        buf = fig2data ( fig )
        w, h, d = buf.shape
        return Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
    
    all_files_set = set(all_files)
    figure_dpi = 0.
    radial_images = {}
    for name1 in all_files_set:
        for name2 in all_files_set:
            try:
                fig = radial_graph_creation(overlap_all[name1,name2][0:3],name1,name2,overlap_all[name1,name2][3],overlap_all[name1,name2][4])
                radial_images[name1,name2] = fig2img(fig)
                figure_dpi = fig.dpi
            except KeyError:
                pass
   
    def bar_graph_creation(pop_tuples):
        fig = plt.figure(figsize=(4,1.5))
        fig.set_facecolor('white')
        ax = fig.add_subplot(111)
        values = {'A':[],'B':[]}
        complementary_tags = {"A":"B","B":"A","M":"M"}
        types = []
        acc = 0
        
        for pt in pop_tuples:
            if len(pt) == 2:
                values[pt[1]].append(pt[0])
                values[complementary_tags[pt[1]]].append(0)
                types.append(pt[1])
                acc += pt[0]
            else:
                values['A'].append(pt[0])
                values['B'].append(pt[1])
                types.append(pt[2]) # 'M'
                acc += pt[0] + pt[1]
        
        x = numpy.array(range(len(values['A'])))
        all_colors = {"A":'#468C54',"B":'#4786A1',"M":'#949E5D','MA':'#468C54','MB':'#4786A1'}
        
        this_colors = []
        for i in range(len(values)):
            if types[i] != 'M':
                this_colors.append(all_colors[types[i]])
            else:
                this_colors.append(all_colors[types[i]+'A'])
        ax.bar(left = x,width = 1.0, height = values['A'],color=this_colors)
        
        this_colors = []
        for i in range(len(values)):
            if types[i] != 'M':
                this_colors.append(all_colors[types[i]])
            else:
                this_colors.append(all_colors[types[i]+'B'])
        ax.bar(left = x, width = 1.0, bottom = values['A'], height = values['B'],color=this_colors)
        
        ax.set_ylabel('%')
        yticks(fontsize=9)
        xticks(fontsize=9)
        fig.text(0.75,0.8,"Acc: %0.2f"%(acc),size= 7)
        return fig
    
    pop_plot_images = {}
    for names in cluster_populations_all:
        fig = bar_graph_creation(cluster_populations_all[names])
        pop_plot_images[names] = fig2img(fig)
        figure_dpi = fig.dpi
    
    margin = 2
    complete_margin = margin*2
    number_of_tiles = len(all_files_set)
    tile_dimensions={'x':int(4*figure_dpi),'y':int(5.5*figure_dpi)+3}
    pop_tile_dimensions = {'x':int(4*figure_dpi),'y':int(1.5*figure_dpi)+3}
    background_dimensions = {'x':int(tile_dimensions['x']+complete_margin)*number_of_tiles,'y':int(tile_dimensions['y']+complete_margin)*number_of_tiles}
    population_offset = int(4*figure_dpi)
    blank_image = Image.new("RGBA", (background_dimensions['x'],background_dimensions['y']))
    grey_tile = Image.new("RGBA", (tile_dimensions['x'],tile_dimensions['y']),color=(220,)*4)
    white_tile = Image.new("RGBA", (tile_dimensions['x'],tile_dimensions['y']),color=(255,)*4)
    grey_pop_tile = Image.new("RGBA", (pop_tile_dimensions['x'],pop_tile_dimensions['y']),color=(220,)*4)
    all_file_names = list(all_files_set)
    for i in range(len(all_file_names)):
        for j in range(len(all_file_names)):
            r_location = (margin+(i*(tile_dimensions['x']+complete_margin)),(j*(tile_dimensions['y']+complete_margin)))
            blank_image.paste(white_tile, r_location)
            
            # Radial graph
            try: 
                blank_image.paste(radial_images[all_file_names[i],all_file_names[j]], r_location)
            except KeyError:
                blank_image.paste(grey_tile, r_location) 
            
            #Population graph
            p_location = (r_location[0],r_location[1]+population_offset)
            try: 
                blank_image.paste(pop_plot_images[all_file_names[i],all_file_names[j]], p_location)
            except KeyError:
                blank_image.paste(grey_pop_tile, p_location) 
            
    blank_image.show()   
    
    for name in all_clustering_sizes:
        print name,numpy.mean(all_clustering_sizes[name]), numpy.std(all_clustering_sizes[name]), all_clustering_sizes[name]
        
    