'''
Created on 13/08/2012

@author: victor
'''
import prody
from pyproclust.clustering.clusterization import Clustering, Cluster
from pyproclust.tools.pdbTools import getPDBStructure



class PCAMetric(object):
    def __init__(self):
        prody.setVerbosity('none')
    
    def evaluate(self,clustering,pdb_struct):
        # Now do calculation for each one of the clusters
        pca_mean_val = 0.;
        for c in clustering.clusters:
            ensemble = prody.Ensemble('pcametric_ensemble')
            ensemble.setCoords( pdb_struct.getCoords())
            ensemble.addCoordset( pdb_struct.getCoordsets()[c.all_elements] )# -> this will work as long as the coordsets are numpy arrays
            ensemble.iterpose()
            pca = prody.PCA('pcametric_pca')
            pca.buildCovariance(ensemble)
            pca.calcModes(n_modes=1)
            pca_mean_val += pca.getEigvals()[0]*c.get_size()
        return  pca_mean_val /clustering.total_number_of_elements
     
if __name__ == "__main__":
    calculator = PCAMetric()
    clustering = Clustering([Cluster(1,range(0,100)),Cluster(100,range(100,200)),Cluster(200,range(200,300))])
    pdb = getPDBStructure("merged_and_extracted_5_0_75.pdb", "resnum < 71")
    calculator.evaluate(clustering, pdb)