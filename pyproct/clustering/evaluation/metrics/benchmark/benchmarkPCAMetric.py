"""
Created on 21/12/2012

@author: victor
"""
import time
import os
import prody
import numpy 
import bz2
from pyRMSD.utils.proteinReading import flattenCoords
from pyproct.clustering.metrics.pcaMetrics import PCAMetric
if __name__ == '__main__':
    """
    Compares Prody and pyProClust implementation.
    """
    
    ######################
    # BENCHMARKING
    ######################
    prody.confProDy(verbosity='none')#setVerbosity('none')
    print "Loading file..."
    t1 = time.time()
    print "\tUncompressing..."
    open("tmp_amber_long.pdb","w").write(bz2.BZ2File("data/amber_long.pdb.tar.bz2").read())
    print "\tLoading..."
    pdb = prody.parsePDB("tmp_amber_long.pdb", subset='calpha')
    not_iterposed_coordsets = numpy.array(pdb.getCoordsets())
    number_of_conformations = not_iterposed_coordsets.shape[0]
    atoms_per_conformation = not_iterposed_coordsets.shape[1]
    os.system("rm tmp_amber_long.pdb")
    print "\tDeleting temporary file"
    t2 = time.time()
    print 'Loading took %0.3f s' % (t2-t1)
    
    ######################
    # PRODY
    ######################
    print "Performing calculations with prody..."
    t1 = time.time()
    ensemble = prody.Ensemble('pca_test_ensemble')
    ensemble.setCoords( pdb.getCoords())
    ensemble.addCoordset(pdb.getCoordsets())
#     prody.setVerbosity('info')
    ensemble.iterpose()
    coordsets = ensemble.getCoordsets()
    pca = prody.PCA('pcametric_pca')
    pca.buildCovariance(ensemble)
    pca.calcModes(n_modes=10)
    t2 = time.time()
    print 'It took %0.3f s' % (t2-t1)
    print "Biggest Eigenvalue: ", pca.getEigvals()
    
    ######################
    # PYPROCLUST
    ######################
    print "Performing calculations with pyProClust..."
    t1 = time.time()
    fcoords = flattenCoords(not_iterposed_coordsets)
    PCAMetric.do_iterative_superposition(fcoords, pdb.getCoordsets().shape[0], pdb.getCoordsets().shape[1])
#     fcoords_copy = fcoords.reshape((number_of_conformations, atoms_per_conformation, 3))
    my_cov_matrix = PCAMetric.create_covariance_matrix(not_iterposed_coordsets)
    biggest_eigenvalue = PCAMetric.calculate_biggest_eigenvalue(my_cov_matrix)
    t2 = time.time()
    print 'It took %0.3f s' % (t2-t1)
    print "Biggest Eigenvalue: ", biggest_eigenvalue
    
