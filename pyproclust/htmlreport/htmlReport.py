'''
Created on 23/10/2012

@author: victor
'''
import numbers
import datetime
#from pyproclust.htmlreport.javascript import page_footer_chunk, page_header_chunk
from pyproclust.tools.scriptTools import create_directory
from PIL import Image
import os

class HTMLReport(object):
    def __init__(self):
        self.report = {'Tries':{
                              'Contents':{ 
                                          'Number of tries': 0,
                                          'All Clusterings': None
                                          },
                              'KO Clusterings': [],
                              'OK Clusterings': []
                              },
                       
                     'Trajectories' : None, 
                     
                     "Matrix Handler": None,
                     
                     'Evaluations': "",
                     
                     'Scores':[],
                     
                     'Best Clustering Selection': [],
                     
                     'Timing':"",
                     
                     "KL" : None,
                     
                     "Image Paths": {
                                     "kl" : "",
                                     "matrix" : "",
                                     "clustering_small":"",
                                     "clustering_big":""
                                     }
                     }
    
    def stringify_analysis_result(self,number):
        if isinstance(number, numbers.Integral):
            return "%d"%number
        elif isinstance(number, float):
            return "%.3f"%number
        else:
            return str(number)
    
    def order_analysis_keys(self,mykeys):
        top_keys = ["Details","Number of clusters","Mean cluster size","Noise level"]
        ordered_keys = []
        
        last_keys = []
        for key in mykeys:
            if "Score" in key:
                last_keys.append(key)
        
        middle_keys = []
        for key in mykeys:
            if not key in last_keys and not key in top_keys:
                middle_keys.append(key)
        
        
        for key in top_keys+middle_keys+last_keys:
            ordered_keys.append(key)
        
        return ordered_keys
        
    def generateHTML(self):
        lines = []
        now = datetime.datetime.now()
        
        more_string = '<a class = "clik_for_details"> (more)</a>'
        
#        "["+now.strftime('%Y-%m-%d %H:%M")+"]
        
        lines.append('<div class="listControl" ><a id="expandList">Expand All</a><a id="collapseList">Collapse All</a></div>')
        
        ##################
        # TABLE HEADER SECTION
        ##################
        lines.append('<div id="runDetailsListContainer" class="listContainer">')
        lines.append('<ul id="runDetailsList" class="expList"> <br>Run details<br> ') 
        ##################
        #TRAJECTORIES SECTION
        ##################
        trajectoryHandler = self.report["Trajectories"]
        lines.append("<li> Trajectory details "+more_string)
        if trajectoryHandler.pdb2 != 0:
            lines.append("<ul> <li>We have compared two pdb files:<br>")
            lines.append("&emsp;&emsp;"+trajectoryHandler.pdb1+" (A) with %d structures"%(trajectoryHandler.pdb1_number_of_conformations)+"<br>")
            lines.append("&emsp;&emsp;"+trajectoryHandler.pdb2+" (B) with %d structures"%(trajectoryHandler.pdb2_number_of_conformations))
        else:
            lines.append("<ul> <li>We have clustered the pdb file:<br>")
            lines.append("&emsp;"+trajectoryHandler.pdb1+" with %d structures"%(trajectoryHandler.pdb1_number_of_conformations))
        lines.append('</li>')
        lines.append(" </ul>")
        lines.append(" <br>")
        lines.append('</li>')
        
        ##################
        # TIMING SECTION
        ##################
        lines.append("<li> Timing "+more_string)
        timing_lines = self.report['Timing'].split("\n")
        lines.append("<ul>")
        for tl in timing_lines:
            if len(tl)>4:
                lines.append("<li> "+tl+"</li>")
        lines.append('</li>')
        lines.append("</ul>")
        lines.append("</div>")
        
        
        lines.append('<div id="summaryListContainer" class="listContainer">')
        lines.append('<ul id="runSummaryList" class = "expList"> <br>Run summary<br> ') 
        ##################
        # FIRST LINE SECTION
        ##################
        lines.append("<li> We have produced "+str(self.report["Tries"]["Contents"]["Number of tries"])+" clusterings "+more_string)
        lines.append("<ul>")
        tags, counter = self.report['Tries']['Contents']['All Clusterings']
        total = 0
        for key in tags:
            lines.append( "<li>"+str(counter[key])+" are of the type "+key+".</li>")
            total += counter[key]
        lines.append("</br>")
        lines.append("</ul>")
        lines.append("</li>")
        
        ##################
        # KO TRIES SECTION
        ##################
        lines.append("<li> "+str(total)+" clusterings were send to the filtering pipeline "+more_string)
        if len(self.report['Tries']['KO Clusterings']) > 0:
            lines.append("<li> Of these, "+str(len(self.report['Tries']['KO Clusterings']))+" were automatically rejected because their number of\
            clusters or the ammount of noise")
            lines.append("<ul>")
            for clustering, reasons in self.report['Tries']['KO Clusterings']:
                lines.append("<li> Rejection details " + more_string)
                lines.append("<ul>")
                lines.append(clustering.details)
                lines.append("Reasons: ")
                lines.append(reasons)
                lines.append("</li>")
            lines.append("</br>")
            lines.append("</ul>")
            lines.append("</li>")
        else:
            lines.append("<ul> No clusterings were rejected because of its number of clusters or its noise</ul><br>")
        lines.append("</li>")
        ##################
        # OK TRIES SECTION
        ##################
        lines.append("<li>"+str(len(self.report['Tries']['OK Clusterings']))+" clusterings were send to the evaluation pipeline "+more_string)
        lines.append("<ul>")
        for clustering in self.report['Tries']['OK Clusterings']:
            lines.append("<li>"+clustering.details+"</li>")
        lines.append("</ul>")
        lines.append("<br>")
        lines.append("</li>")
        
        ##################
        # EVAL SECTION
        ##################
        # Do a clustering->analysis_dic dictionary
        results_pack = self.report["Evaluations"]
        all_analysis_dic = {}
        for clustering,analysis_dic in results_pack:
            all_analysis_dic[clustering] = analysis_dic
        (best_score, best_clustering) = self.report["Best Clustering Selection"]
        # Then index by scoring
        all_scores = self.report['Scores']
        lines.append("<li> We used "+str(len(self.report['Scores']))+" different criteria to do evaluations "+more_string)
        i = 0
        lines.append("<ul>")
        for value_map, scores in all_scores:
            crit_id = "Criteria "+str(i)
            lines.append("<li>"+crit_id+": "+str(value_map)+"</li>")
            for norm_score, clustering in  scores:
                all_analysis_dic[clustering][crit_id+" Score"] = norm_score
            i=  i+1
        lines.append("</ul>")
        lines.append("<br>")
        lines.append("</li>")
        
        ##################
        # EVAL TABLE SECTION
        ##################
        lines.append('<li> Clustering evaluation '+more_string)
        lines.append("<ul>")
        lines.append("<li><br><br><br><br>")
        ordered_analysis_keys =  self.order_analysis_keys(all_analysis_dic[all_analysis_dic.keys()[0]].keys())
        lines.append('<table id = "summary_table" >\n')
        lines.append('<thead>')
        lines.append('<tr>')
        for key in ordered_analysis_keys:
            if key == "Details":
                lines.append('<th class = "details"><b>'+str(key)+"</b></td>\n")
            else:
                lines.append('<th class="vertical"><b>'+str(key)+"</b></td>\n")
        lines.append("</tr>\n")
        lines.append('</thead>')
        lines.append('<tbody>')
        for clustering in all_analysis_dic.keys():
            if clustering == best_clustering:
                detail_class = 'class ="best_clustering_details"'
                normal_class =  'class ="best_clustering"'
            else:
                detail_class = 'class ="details"'
                normal_class =  ''
            lines.append("<tr "+normal_class+">\n")
            
            for key in ordered_analysis_keys:
                if key == "Details":
                    lines.append('<td '+detail_class+'><b>'+self.stringify_analysis_result(all_analysis_dic[clustering][key])+"</b></td>\n")
                else:
                    if key == "Noise level":
                        lines.append("<td "+normal_class+">"+self.stringify_analysis_result(all_analysis_dic[clustering][key]*100)+"%</td>\n")
                    else:
                        lines.append("<td "+normal_class+">"+self.stringify_analysis_result(all_analysis_dic[clustering][key])+"</td>\n")
            lines.append("</tr>\n")
        lines.append('</tbody>')
        lines.append('''<tfoot>
        <tr>
        <td><div class = "green_square"></div><div class="green_square_text">Best clustering</div></td>
        </tr>
        </tfoot>''')
        lines.append("</table>\n")
        lines.append("<li>")
        lines.append("</ul>")
        lines.append('</li>')    
        lines.append('</ul>')
        lines.append("</div>")
        
        ##############################
        # GRAPHICAL RESULTS SECTION
        ##############################
        lines.append('<div id="resultsListContainer" class="listContainer">')
        lines.append('<ul id="resultsList" class ="expList"> <br>Results<br> ')
        lines.append('<li>RMSD Matrix '+more_string+'<ul><li><img class="thumbnail" src="img/thumbnails/matrix_plot_thumbnail.jpg"></li>')
        matrixHandler = self.report["Matrix Handler"]
        lines.append("<li><table>")
        lines.append("<tr><td>Minimum:</td><td>%.5f</td><tr>"%(matrixHandler.min_dist ))
        lines.append("<tr><td>Maximum:</td><td>%.5f</td><tr>"%(matrixHandler.max_dist ))
        lines.append("<tr><td>Mean:</td><td>%.5f</td><tr>"%(matrixHandler.mean_dist ))
        lines.append("<tr><td>Std. Dev.:</td><td>%.5f</td><tr>"%(matrixHandler.std))
        lines.append("<tr><td>Skewness:</td><td>%.5f</td><tr>"%(matrixHandler.skew))
        lines.append("<tr><td>Kurtosis:</td><td>%.5f</td><tr>"%(matrixHandler.kurtosis))
        lines.append("</table></li></ul></li>")
        
        lines.append('<li>RMSD probability distributions'+more_string+'<ul><li>')
        klDiv = self.report["KL"]
        lines.append('<a href="http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence">Kullback - Leibler divergence </a>:<br>')
        lines.append(' &emsp;&emsp;A - B: %.3f<br>'%klDiv.kl1)
        lines.append(' &emsp;&emsp;B - A: %.3f<br>'%klDiv.kl2)
        lines.append('</li><li><img class="thumbnail" src="img/thumbnails/rmsd_distrib_thumbnail.jpg"></li></ul></li>')
        
        lines.append('<li>Clustering Results (small) '+more_string+'<ul><li><img class="thumbnail" src="img/thumbnails/analysis_plot_small_thumbnail.jpg"></li></ul></li>')
        
        lines.append('<li>Clustering Results (big) '+more_string+'<ul><li><img class="thumbnail" src="img/thumbnails/analysis_plot_big_thumbnail.jpg"></li></ul></li>')
        lines.append('</ul>')
        lines.append("</div>")
        
        final_html = ""#page_header_chunk
        for l in lines:
            final_html += l+"\n"
        
        #final_html += page_footer_chunk 
        
        return final_html
    
    def create_thumbnails(self, dest_path):
        # Create img path there
        create_directory(dest_path+"/img/thumbnails")
        thumb_size = (256,256)
        for key in self.report["Image Paths"]:
            file_w_all, ext = os.path.splitext(self.report["Image Paths"][key])
            file_name = file_w_all.split("/")[-1]
            im = Image.open(self.report["Image Paths"][key])
            im.thumbnail(thumb_size, Image.ANTIALIAS)
            im.save(dest_path+"/img/thumbnails/"+file_name + "_thumbnail.jpg"  , "JPEG")
        