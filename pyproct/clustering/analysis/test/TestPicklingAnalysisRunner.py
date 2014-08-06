"""
Created on 15/06/2012

@author: victor
"""
import unittest
import pickle
import tempfile
from pyproct.clustering.analysis.picklingAnalysisRunner import PicklingAnalysisRunner,\
    run_all_analysis_for_a_clustering
from pyproct.clustering.analysis.analysis import Analysis

class ClusteringMock():
    def __init__(self, name, number_of_clusters, number_of_elements):
        self.number_of_clusters = number_of_clusters
        self.total_number_of_elements = number_of_elements
        self.clusters = range(number_of_clusters)
        self.details = "ClusteringMock"
        self.name = name
        
    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.name
        
class SchedulerMock():
    def __init__(self):
        self.queue = []
        
    def add_process(self, process_name, description, run_all_analysis_for_a_clustering, func_kwargs):
        self.queue.append((process_name, description, run_all_analysis_for_a_clustering, func_kwargs))
        
class PopulatorMock():
    def __init__(self):
        pass
    
    def get_analysis_list(self):
        return []

class TestPicklingAnalysisRunner(unittest.TestCase):
    
    def test_run_analysis_for(self):
        scheduler = SchedulerMock()
        prunner = PicklingAnalysisRunner(scheduler, "", {}, PopulatorMock())
        
        prunner.run_analysis_for_a_clustering("clustering1",
                                              {
                                               "clustering1":{
                                                              "type":"algorithm2",
                                                              "clustering":ClusteringMock("clustering1", 100, 10000),
                                                              "parameters":{}
                                                              },
                                               "clustering2":{
                                                              "type":"algorithm1",
                                                              "clustering":ClusteringMock("clustering2", 50, 5000),
                                                              "parameters":{}
                                                              }
                                               })
        
        # Assert that the process has been scheduled
        self.assertEqual(prunner.current_analysis, 1)
        self.assertEqual(len(scheduler.queue), 1)
        process_name, description, unused1, unused2 = scheduler.queue.pop()
        self.assertEqual(process_name, 'Evaluation of clustering1')
        self.assertEqual(description, 'Evaluation of algorithm2 clustering (clustering1) with parameters: {}')
        
        scheduler = SchedulerMock()
        prunner = PicklingAnalysisRunner(scheduler, "", {}, PopulatorMock())
        prunner.run_analysis_for_all_clusterings({
                                               "clustering1":{
                                                              "type":"algorithm2",
                                                              "clustering":ClusteringMock("clustering1", 100, 10000),
                                                              "parameters":{}
                                                              },
                                               "clustering2":{
                                                              "type":"algorithm1",
                                                              "clustering":ClusteringMock("clustering2", 50, 5000),
                                                              "parameters":{}
                                                              }
                                               })
        
        self.assertEqual(prunner.current_analysis, 2)
        self.assertEqual(len(scheduler.queue), 2)
        process_name1, description1, unused1, unused2 = scheduler.queue[0]
        self.assertEqual(process_name1, 'Evaluation of clustering1')
        self.assertEqual(description1, 'Evaluation of algorithm2 clustering (clustering1) with parameters: {}')
        process_name2, description2, unused1, unused2 = scheduler.queue[1]
        self.assertEqual(process_name2, 'Evaluation of clustering2')
        self.assertEqual(description2, 'Evaluation of algorithm1 clustering (clustering2) with parameters: {}')
        
         
    def test_run_and_recover_data(self):
        
        def analysis(clustering, adjective):
            return "Analysis %s performed for %s"%(adjective, str(clustering))
        
        analysis_queue = []
        analysis_queue.append( Analysis("First analysis", analysis, "First"))
        analysis_queue.append( Analysis("Second analysis", analysis, "Second"))
        
        clustering_ids = ['Clustering 1','Clustering 3','Clustering 2', 'Clustering 4']
         
        scheduler = SchedulerMock()
        prunner = PicklingAnalysisRunner(scheduler, "", {}, PopulatorMock())

        for clustering_id in clustering_ids:
            tmp_file_handler_descriptor , path = tempfile.mkstemp() 
            prunner.evaluation_data.append(path)
            run_all_analysis_for_a_clustering(clustering_id, 
                                              ClusteringMock(clustering_id, 100, 10000), 
                                              tmp_file_handler_descriptor, 
                                              analysis_queue)

        clustering_info = {
                            'Clustering 1':{
                            },
                            'Clustering 3':{                                            
                            },
                            'Clustering 2':{                                            
                            },
                            'Clustering 4':{
                            }
        }
        
        prunner.recover_evaluation_data(clustering_info)
        
        expected = {
                    'Clustering 1': {
                                     'evaluation': {
                                                    'First analysis': 'Analysis First performed for Clustering 1',
                                                    'Second analysis': 'Analysis Second performed for Clustering 1'
                                                    }
                    },
                    'Clustering 2': {
                                     'evaluation': {
                                                    'First analysis': 'Analysis First performed for Clustering 2',
                                                    'Second analysis': 'Analysis Second performed for Clustering 2' 
                                                    }
                    }, 
                    'Clustering 3': {
                                     'evaluation': {
                                                    'First analysis': 'Analysis First performed for Clustering 3',
                                                    'Second analysis': 'Analysis Second performed for Clustering 3'
                                                    }
                    }, 
                    'Clustering 4': {
                                     'evaluation': {
                                                    'First analysis': 'Analysis First performed for Clustering 4',
                                                    'Second analysis': 'Analysis Second performed for Clustering 4'
                                                    }
                    } 
        }
        
        self.assertDictEqual(expected,clustering_info)
          
    def test_temporary_file_check(self):
        myobject = ["Clustering 1", "Clustering 3", "Clustering 2","Clustering 4"]
        tmp_file_handler = tempfile.TemporaryFile()
        pickle.dump(myobject, tmp_file_handler)
        tmp_file_handler.seek(0)
        self.assertEqual(myobject, pickle.load(tmp_file_handler))
         
#     def test_repack_results(self):
#         prunner = PicklingParallelAnalysisRunner(100)
#         ordered_clusterings = ["Clustering 1", "Clustering 3", "Clustering 2","Clustering 4"]
#         all_results = {"number":[1,2,3,4],"color":["red","green","blue","gray"]}
#         expected_repack = [('Clustering 1', {'color': 'red', 'number': 1}), ('Clustering 3', {'color': 'green', 'number': 2}), ('Clustering 2', {'color': 'blue', 'number': 3}), ('Clustering 4', {'color': 'gray', 'number': 4})]
#         self.assertItemsEqual(expected_repack, prunner.repack_results(ordered_clusterings, all_results))
#          
#     def test_gen_final_string_and_normalize(self):
#         prunner = PicklingParallelAnalysisRunner(100)
#         all_results = {"number":[1,2,3,4],"color":["red","green","blue","gray"]}
#         expected_result_string = "color\tred\tgreen\tblue\tgray\t\nnumber\t1\t2\t3\t4\t\n"
#         self.assertEqual(expected_result_string, prunner.gen_final_string_and_normalize(all_results))
#         self.assertItemsEqual(all_results["number"],[ 0.25, 0.5, 0.75, 1.])
#         self.assertItemsEqual(all_results["color"],["red","green","blue","gray"])
# 
#     def test_gen_results_string(self):
#         prunner = PicklingParallelAnalysisRunner(100)
#         self.assertEqual("number\t1\t2\t3\t4\t", prunner.gen_results_string( "number", [1,2,3,4]))
#         self.assertEqual("color\tred\tgreen\tblue\t", prunner.gen_results_string( "color", ["red","green","blue"]))
#     
#     def test_numerical_results_normalization(self):
#         arr = [1,2,3,4,5,6,7,8,9,10]
#         expected_result = [ 0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1. ]
#         prunner = PicklingParallelAnalysisRunner(100)
#         numpy.testing.assert_array_equal(expected_result, prunner.numerical_results_normalization(arr))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_generate_report']
    unittest.main()