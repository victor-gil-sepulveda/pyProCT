'''
Created on 29/04/2013

@author: victor
'''
import json
from pyproct.clustering.cluster import Cluster
from pyproct.clustering.clustering import Clustering

#http://stackoverflow.com/questions/4821940/how-to-make-simplejson-serializable-class
class SerializerRegistry(object):
    def __init__(self):
        self._classes = {}
    
    def add(self, cls):
        self._classes[cls.__module__, cls.__name__] = cls
        return cls
    
    def object_hook(self, dct):
        module, cls_name = dct.pop('__type__', (None, None))
        if cls_name is not None:
            return self._classes[module, cls_name].from_dic(dct)
        else:
            return dct
        
    def default(self, obj):
        return obj.to_dic()

class ClusteringResultsGatherer(object):
    def __init__(self):
        pass
    
    def gather(self, timer_handler, trajectory_handler, workspace_handler, clustering_results, files):
        results = {}
        results["timing"] = timer_handler.get_elapsed()
        results["trajectories"] = trajectory_handler.pdbs
        if(clustering_results is not None):
            results["best_clustering"] = clustering_results[0]
            results["selected"] = clustering_results[1]
            results["not_selected"] = clustering_results[2]
            results["scores"] = clustering_results[3]
        results["files"] = files
        results["workspace"] = workspace_handler.data

        serializer = SerializerRegistry()
        serializer.add(Clustering)
        serializer.add(Cluster)
        return json.dumps(results, 
                          sort_keys=False, 
                          indent=4, 
                          separators=(',', ': '),
                          default=serializer.default)