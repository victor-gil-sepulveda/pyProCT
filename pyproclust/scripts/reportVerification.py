import sys
import pickle
from pyproclust.clustering.selection.bestClusteringSelector import BestClusteringSelector
import numpy
import matplotlib.pyplot as plt
import scipy.stats

###################################
#
# Workbench for parameter estimation.
#
###################################

if __name__ == '__main__':
    report_file = sys.argv[1]
    reports_file_handler = open(report_file,"r")
    results_pack = pickle.load(reports_file_handler)
    reports_file_handler.close()
    
    EXPECTED_NUMBER_OF_CLUSTERS = int(sys.argv[2]) # 3
    ELEMENTS_PER_CLUSTER = int(sys.argv[3]) # 100
    
    NUM_CLUSTERS_CUTOFF = 20
    
    
    cluster_score_value_map = [{"CythonSilhouette":(0.9,">"),"PCAanalysis":(1.,"<"),"NormNCut":(0.3,">")},\
                               {"CythonSilhouette":(0.8,">"),"PCAanalysis":(1.,"<"),"CythonMirrorCohesion":(0.2,"<"),"NormNCut":(0.1,">")}]
    
    #########################
    ## Adding a new analysis
    #########################
    def compare_with_golden(clustering):
        """
        Calculate % of success (we know that in our examples there are EXPECTED_NUMBER_OF_CLUSTERS [N] 
        clusters of ELEMENTS_PER_CLUSTER elements 
        which would be like this [0,99], [100,199], ... [N-100,N-1]
        """
        clusters = clustering.clusters
        if len(clusters) != EXPECTED_NUMBER_OF_CLUSTERS:
            return 0.
        ranges = []
        for i in range(EXPECTED_NUMBER_OF_CLUSTERS):
            ranges.append(set(range(i*ELEMENTS_PER_CLUSTER,(i+1)*ELEMENTS_PER_CLUSTER)))
        # Do some mapping
        cluster_map = {} 
        for cluster in clusters:
            elems = cluster.all_elements
            # To which range does this belong?
            elems_set = set(elems)
            card_intersection = []
            for j in range(EXPECTED_NUMBER_OF_CLUSTERS):
                card_intersection.append((len(elems_set.intersection(ranges[j])),j))
            card_intersection.sort(reverse=True)
            similarity_score = card_intersection[0][0] #@UnusedVariable
            range_pert = card_intersection[0][1]
            cluster_map[range_pert] = similarity_score
        # If we are not able to do a mapping...
        if len(cluster_map.keys()) != EXPECTED_NUMBER_OF_CLUSTERS:
            return 0.
        # If we can, let's see how much misses we did (percentually)
        scores = []
        for k in cluster_map:
            scores.append(cluster_map[k])
        return numpy.mean(scores) / float(ELEMENTS_PER_CLUSTER)
    
    
    for (clustering, analysis_dic) in results_pack:
        analysis_dic["NormNCut"] = analysis_dic["NCut"]/len(clustering.clusters)
        analysis_dic["GoldenComp"] =  compare_with_golden(clustering)
        print clustering.details, compare_with_golden(clustering), len(clustering.clusters)
#        if "Spectral Clustering (k = 9)" in clustering.details:
#            for c in clustering.clusters:
#                print sorted(c.all_elements) 
    
    bcSel = BestClusteringSelector(cluster_score_value_map)
    
    best_score, best_clustering =  bcSel.chooseBestClustering(results_pack) 
    
    print "Winner is ", best_clustering.details, " with ", best_score, "(",len(best_clustering.clusters)," clusters) and correctness: ",compare_with_golden(best_clustering) 
    

    #allowed_analysis = ["NormNCut","CythonSilhouette","NCut","MixMaxCut","RatioCut","CythonMirrorCohesion","PCAanalysis","Number of clusters","GoldenComp","CythonMinimumMeanSeparation"]
    allowed_analysis = ["NormNCut","CythonSilhouette","NCut","MixMaxCut","CythonMirrorCohesion","PCAanalysis","GoldenComp","CythonMinimumMeanSeparation"]
    all_analysis_dic = {}
    for (clustering, analysis_dic) in results_pack:
        if len(clustering.clusters) < NUM_CLUSTERS_CUTOFF:
            for analysis in analysis_dic:
                if not analysis in all_analysis_dic:
                    all_analysis_dic[analysis] = []
                if analysis == "Number of clusters":
                    all_analysis_dic[analysis].append(len(clustering.clusters))
                else:
                    all_analysis_dic[analysis].append(analysis_dic[analysis])
    

#    for analysis in all_analysis_dic:
#        if analysis in allowed_analysis:
#            print analysis, 
#            for value in all_analysis_dic[analysis]:
#                try:
#                    print "\t%.4f"%value,
#                except TypeError:
#                    print "\t"+value,
#            print
    
    
            
    
    # Plot with matplotlib to see correlation
    #allowed_analysis = ["PCAanalysis","Number of clusters","GoldenComp"]
    
    fig = plt.figure()
    host = fig.add_subplot(111)
    par1 = host.twinx()   
    for analysis in all_analysis_dic:
        p = None
        if analysis in allowed_analysis and analysis != "Number of clusters":
            host.plot(all_analysis_dic[analysis],label = analysis)
        if analysis == "Number of clusters":
            par1.plot(all_analysis_dic[analysis],label = analysis)
    
    handles, labels = host.get_legend_handles_labels()
#    handles2, labels2 = par1.get_legend_handles_labels()
#    handles.extend(handles2)
#    labels.extend(labels2)
    plt.legend(handles,labels)
#    plt.show()
    
    correlations = []
    correlation_pairs = []
    all_analysis_dic_keys = all_analysis_dic.keys()
#    for i in range(len(all_analysis_dic)-1):
#        analysis1 = all_analysis_dic_keys[i] 
#        for j in range(i+1,len(all_analysis_dic)):
#            analysis2 = all_analysis_dic_keys[j]
#            #if analysis1 in allowed_analysis and analysis2 in allowed_analysis:
#            if (analysis1 == "GoldenComp") or (analysis2 == "GoldenComp"):
#                correlations.append((scipy.stats.spearmanr(all_analysis_dic[analysis1],all_analysis_dic[analysis2]), analysis1,analysis2)) 
    for analysis in allowed_analysis:
        if analysis != "GoldenComp":
            correlations.append((scipy.stats.spearmanr(all_analysis_dic[analysis],all_analysis_dic["GoldenComp"]), analysis,"GoldenComp"))
    
    correlations.sort()
    for c in correlations:
        print c
