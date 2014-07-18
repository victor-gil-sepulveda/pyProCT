pyProCT
=======

pyProCT is an open source cluster analysis software especially adapted
for jobs related with structural proteomics. Its approach allows users
to define a clustering goal (clustering hypothesis) based on their
domain knowledge. This hypothesis will guide the software in order to
find the best algorithm and parameters (including the number of
clusters) to obtain the result that better fulfills their expectatives.
In this way users do not need to use cluster analysis algorithms as a
black box, which will (hopefully) improve their results. pyProCT not
only generates a resulting clustering, it also implements some use cases
like the extraction of representatives or trajectory redundance
elimination.

`Table of Contents <http://doctoc.herokuapp.com/>`_

-  `pyProCT <#user-content-pyproct>`_

   -  `Installation <#user-content-installation>`_
   -  `Using pyProCT as a standalone
      program <#user-content-using-pyproct-as-a-standalone-program>`_

      -  `Global <#user-content-global>`_
      -  `Data <#user-content-data>`_
      -  `Clustering <#user-content-clustering>`_

         -  `generation <#user-content-generation>`_
         -  `algorithms <#user-content-algorithms>`_
         -  `evaluation <#user-content-evaluation>`_

      -  `Postprocessing <#user-content-postprocessing>`_
      -  `Checking the script <#user-content-checking-the-script>`_
      -  `Learn more <#user-content-learn-more>`_

   -  `Using pyProCT as part of other
      programs <#user-content-using-pyproct-as-part-of-other-programs>`_

      -  `Using it as a separate program from other Python
         script <#user-content-using-it-as-a-separate-program-from-other-python-script>`_

   -  `Parallel execution <#user-content-parallel-execution>`_

-  `Documentation <#user-content-documentation>`_
-  `TODO <#user-content-todo>`_

Installation
------------

pyProCT is quite easy to install using *pip*. Just write:

``Shell > sudo pip install pyProCT`` And *pip* will take care of all the
dependencies (shown below).

 It is recommended to install Numpy and Scipy before starting the
installation using your OS software manager. You can try to download and
install them
`manually <http://docs.scipy.org/doc/numpy/user/install.html>`_ if you
dare.

 mpi4py is pyProCT's last dependency. It can give problems when
installing it in OS such as SUSE. If the installation of this last
package is not succesful, pyProCT can still work in Serial and Parallel
(using *multiprocessing*) modes.

Using pyProCT as a standalone program
-------------------------------------

The preferred way to use pyProCT is through a JSON "script" that
describes the clustering task. It can be executed using the following
line in your shell:

::

    > python -m pyproct.main script.json

The JSON script has 4 main parts, each one dealing with a different
aspect of the clustering pipeline. This sections are: \* *"global"*:
Handles workspace and scheduler parameterization. \* *"data"*: Handles
distance matrix parameterization. \* *"clustering"*: Handles algorithms
and evaluation parameterization. \* *"preprocessing"*: Handles what to
do with the clustering we have calculated.

::

    {
        "global":{},
        "data":{},
        "clustering":{},
        "postprocessing":{}
    }

Global
~~~~~~

``JSON {     "control": {         "scheduler_type": "Process/Parallel",         "number_of_processes": 4     },     "workspace": {          "tmp": "tmp",          "matrix": "matrix",          "clusterings": "clusterings",          "results": "results",          "base": "/home/john/ClusteringProject"     } }``
This is an example of *"global"* section. It describes the work
environment (workspace) and the type of scheduler that will be built.
Defining the subfolders of the wokspace is not mandatory, however it may
be convenient in some scenarios (for instance, in serial multiple
clustering projects, sharing the *tmp* folder would lower the disk usage
as at each step it will be overwritten).

This is a valid global section using a serial scheduler and default
names for workspace inner folders:
``JSON {     "control": {         "scheduler_type": "Serial"     },     "workspace": {          "base": "/home/john/ClusteringProject"     } }``
pyProCT allows the use of 3 different schedulers that help to improve
the overall performance of the software by parallelizing some parts of
the code. The available schedulers are "Serial", "Process/Parallel"
(uses Python's
`multiprocessing <https://docs.python.org/2/library/multiprocessing.html>`_)
and "MPI/Parallel" (uses MPI through the module
`mpi4py <http://mpi4py.scipy.org/>`_).

Data
~~~~

The *"data"* section defines how pyProCT must build the distance matrix
that will be used by the compression algorithms. Currently pyProCT
offers up to three options to build that matrix: "load", "rmsd" and
"distance" - *"rmsd"*: Calculates a all vs all rmsd matrix using any of
the
`pyRMSD <https://github.com/victor-gil-sepulveda/pyRMSD#collective-operations>`_
calculators available. It can calculate the RMSD of the fitted region
(defined by `Prody <http://prody.csb.pitt.edu/>`_ compatible selection
string in *fit\_selection*) or one can use one selection to superimpose
and another to calculate the rmsd (*calc\_selection*). There are two
extra *parameters* that must be considered when building an RMSD matrix.
- *"type"*: This property can have two values: *"COORDINATES"* or
*"DIHEDRALS"*. If *DIHEDRALS* is chosen, each element *(i,j)* of the
distance matrix will be the RMSD of the arrays containing the phi-psi
dihedral angle series of conformation *i* and *j*. - *"chain\_map"*: If
set to *true* pyProCT will try to reorder the chains of the biomolecule
in order to minimize the global RMSD value. This means that it will
correctly calculate the RMSD even if chain coordinates were permuted in
some way. The price to pay is an increase of the calculation time
(directly proportional to the number of chains or the number of chains
having the same length). - *"distance"*: After superimposing the
selected region it calculates the all vs all distances of the
geometrical center of the region of interest (*body\_selection*). -
*"load"*: Loads a precalculated matrix.

JSON chunk needed to generate an RMSD matrix from two trajectories:
``JSON {     "type": "pdb_ensemble",     "files": [         "A.pdb",         "B.pdb"     ],     "matrix": {         "method": "rmsd",         "parameters": {             "calculator_type": "QCP_OMP_CALCULATOR",             "fit_selection": "backbone",         },         "image": {             "filename": "matrix_plot"         },         "filename":"matrix"     } }``
JSON chunk to generate a dihedral angles RMSD matrix from one
trajectories:
``JSON {     "type": "pdb_ensemble",     "files": [         "A.pdb"     ],     "matrix": {         "method": "rmsd",         "parameters": {             "type":"DIHEDRAL"         },         "image": {             "filename": "matrix_plot"         },         "filename":"matrix"     } }``

The matrix can be stored if the *filename* property is defined. The
matrix can also be stored as an image if the *image* property is
defined.

pyProCT can currently load *pdb* and *dcd* files. The details to load
the files must be written into the array under the "files" keyword.
There are many ways of telling pyProCT the files that have to be load
and can be combined in any way you like:

1. Using a list of file paths. If the file extension is ".txt" or
   ".list" it will be treated as a pdb list file. Each line of such
   files will be a pdb path or a pdb path and a selection string,
   separated by comma.

``A.pdb, name CA B.pdb C.pdb, name CA ...`` 2. Using a list of file
objects: ``JSON {     "file": ... ,     "base_selection": ... }`` Where
*base\_selection* is a `Prody <http://prody.csb.pitt.edu/>`_ compatible
selection string. Loading files this way can help in cases where not all
files have structure with the same number of atoms: *base\_selection*
should define the common region between them (if a 1 to 1 map does not
exist, the RMSD calculation will be wrong).

3. Only for *dcd* files:
   ``JSON {     "file": ...,     "atoms_file": ...,     "base_selection": ... }``
   Where *atoms\_file* is a *pdb* file with at least one frame that
   holds the atomic information needed by the *dcd* file.

**Note*: data.type is currently unsused*

Clustering
~~~~~~~~~~

The *clustering* section specifies how the clustering exploration will
be done. It is divided in 3 other subsections:

``JSON {     "generation": {         "method": "generate"     },     "algorithms": {         ...     },     "evaluation": {         ...     } }``
#### Generation Defines how to generate the clustering (*"load"* or
*"generate"*). if *"load"* is chosen, this section will also contain the
clustering that may be used in the *"clusters"* property. Ex.:

::

    {
        "clustering": {
            "generation": {
                "method" : "load",
                "clusters": [
                        {
                            "prototype " : 16,
                            "id": "cluster_00",
                            "elements" : "9, 14:20"
                        },
                        {
                            "prototype": 7,
                            "id": "cluster_01",
                            "elements": "0:8, 10:14, 21"
                        }
                ]
            }
    }

Algorithms
^^^^^^^^^^

If clustering.generation.method equals "generate", this section defines
the algorithms that will be used as well as their parameters (if
necessary). The currently available algorithms are : *"kmedoids"*,
*"hierarchical"*, *"dbscan"*, *"gromos"*, *"spectral"* and *"random"*.
Each algorithm can store its list of parameters, however the preferred
way to work with pyProCT is to let it automatically generate them.
Almost all algorithms accept the property *max*, that defines the
maximum amount of parameter collections that will be generated for that
algorithm. Ex.

::

    {
        "kmedoids": {
            "seeding_type": "RANDOM",
            "max": 50,
            "tries": 5
        },
        "hierarchical": {

        },
        "dbscan": {
            "max": 50
        },
        "gromos": {
            "max": 50
        },
        "spectral": {
            "max": 50,
            "force_sparse":true,
        }
    }

Algorithm parameters can be explicitly written, however it is not
recommended:

::

    {
        "kmedoids": {
            "seeding_type": "RANDOM",
            "max": 50,
            "tries": 5,
            "parameters":[{"k":4},{"k":5},{"k":6}]
        }
    }

Evaluation
^^^^^^^^^^

This section holds the *Clustering Hypothesis*, the core of pyProCT.
Here the user can define how the expected clustering will be. First the
user must set the expected number of clusters range. Also, an estimation
of the dataset noise and the cluster minimum size (the minimum number of
elements a cluster must have to not be considered noise) will complete
the quantitative definition of the target result.

Ex.
``JSON {     "maximum_noise": 15,     "minimum_cluster_size": 50,     "maximum_clusters": 200,     "minimum_clusters": 6,     "query_types": [ ... ],     "evaluation_criteria": {         ...     } }``

The second part of the *Clustering Hypothesis* tries to characterize the
clustering internal traits in a more qualitative way. Concepts like
cluster "Compactness" or "Separation" can be used here to define the
expected clustering. To this end users must write their expectations in
form of *criteria*. This criteria are, in general, linear combinations
of Internal Clustering Validation Indices (ICVs). The best clustering
will be the one that gets the best score in any of these *criteria*. See
`this
document <https://dl.dropboxusercontent.com/u/58918851/icv_info.pdf>`_
to get more insight about the different implemented criteria and their
meaning.

Additionally users may choose to ask pyProCT about the results of this
ICVs and other evaluation functions(e.g. the average cluster size) by
adding them to the *queries* array.

::

    {
            ...
        "query_types": [
            "NumClusters",
            "NoiseLevel",
            "MeanClusterSize"
        ],
        "evaluation_criteria": {
            "criteria_0": {
                "Silhouette": {
                    "action": ">",
                    "weight": 1
                }
            }
        }
    }

Postprocessing
~~~~~~~~~~~~~~

Getting a good quality clustering is not enough, we would like to use
them to extract useful information. pyProCT implements some use cases
that may help users to extract this information.

\`\`\`JSON { "rmsf":{},

::

    "centers_and_trace":{},

    "representatives":{
        "keep_remarks": [true/false],
        "keep_frame_number": [true/false]
    },

    "pdb_clusters":{
        "keep_remarks": [true/false],
        "keep_frame_number": [true/false]
    },

    "compression":{
        "final_number_of_frames": INT,
        "file": STRING,
        "type":[‘RANDOM’,’KMEDOIDS’]
    },

    "cluster_stats":{
        "file": STRING
    },


    "conformational_space_comparison":{},

    "kullback_liebler":{}

} \`\`\` - *"rmsf"* : Calculates the global and per-cluster (and
per-residue) root mean square fluctuation (to be visualized using the
`GUI <https://github.com/victor-gil-sepulveda/pyProCT-GUI>`_).

-  *"centers\_and\_trace"* : Calculates all geometrical centers of the
   calculation selection of the system (to be visualized using the
   `GUI <https://github.com/victor-gil-sepulveda/pyProCT-GUI>`_).

-  *"representatives"* : Extracts all the representatives of the
   clusters in the same pdb.
    Parameters:

   -  *"keep\_remarks"*: If true every stored model will be written
      along with its original remarks header. Default: false.
   -  *"keep\_frame\_number"*: If true, the model number of any stored
      conformation will be the original pdb one. Default: false.

-  *"pdb\_clusters"* : Extracts all clusters in separate pdbs.
    Parameters:

   -  *"keep\_remarks"*: If true every stored model will be written
      along with its original remarks header. Default: false.
   -  *"keep\_frame\_number"*: If true, the model number of any stored
      conformation will be the original pdb one. Default: false.

-  *"compression"* : Reduces the redundancy of the trajectory using the
   resulting clustering.
    Parameters:

   -  *"file"*: The name of the output file without extension. Default
      "compressed"(.pdb)
   -  *"final\_number\_of\_frames"*: The expected (minimum) number of
      frames of the compressed file.
   -  *"type"*: The method used to get samples from each cluster.
      Options:

      -  "RANDOM": Gets a random sample of the elements of each cluster.
      -  "KMEDOIDS": Applies the k-medoids algorithm to the elements of
         a cluster and stores the representatives.
         Default: "KMEDOIDS".

-  *"cluster\_stats"*: Generates a human readable file with the
   distances between cluster centers and their diameters.
    Parameters:

   -  *"file"*: The name of the output file without extension (will be
      sotred into the results folder). Default:
      "per\_cluster\_stats"(.csv).

-  *"conformational\_space\_comparison"* : Work in progress.

-  *"kullback\_liebler"* : Work in progress.

Script validation
~~~~~~~~~~~~~~~~~

As the control script is indeed holding a JSON object, any JSON
validator can be used to discover the errors in case of script loading
problems. A good example of such validators is
`JSONLint <http://jsonlint.com/>`_. pyProCT scripts accept javascript
comments ( // and /\* \*/)

 Using pyProCT as part of other programs
----------------------------------------

-  Using algorithms
-  Clustering from label lists
-  Using ICVs with custom clusterings
-  Performing the whole protocol
   ``Python Driver(Observer()).run(parameters)``

The necessary documentation to use pyProCT classes is written inside the
code. It has been extracted
`here <pyproct/docs/_build/html/index.html>`_ and
`here <pyproct/docs/doxyxml/html/index.html>`_. We are currently trying
to improve this documentation with better explanations and examples.

See `this
file <https://github.com/victor-gil-sepulveda/pyProCT/blob/master/validation/bidimensional/validation_main.py>`_.

Using it as a separate program from other Python script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-   Loading results
-   Generating scripts programatically

See `this project <https://github.com/victor-gil-sepulveda/PhD-GPCR>`_
for some examples.

Parallel execution
------------------

To execute pyProCT in parallel you just need to issue this line:

::

    > mpirun -np NumberOfProcesses -m pyproct.main --mpi script.json

 When running pyProCT using MPI you will need to use the *MPI/Parallel*
Scheduler or it will just execute several independent serial runs.

 Remember that you need to use the same libraries and versions to build
mpi4py and mpirun, otherwise you won't be able to execute it.

Documentation
=============

We are still experimenting to see which documentation generator fits
better with us. Currently we have two versions of the documentations:
one using `Sphinx <http://sphinx-doc.org/>`_ and the other using
`Doxygen <http://www.stack.nl/~dimitri/doxygen/>`_+`doxpy <http://code.foosel.org/doxypy>`_.
See them `here <pyproct/docs/_build/html/index.html>`_ and
`here <pyproct/docs/doxyxml/html/index.html>`_. We will possibly publish
it in a cloud solution like
`readthedocs.org <https://readthedocs.org/>`_

Learn more
~~~~~~~~~~

A more detailed explanation of the script contents can be found
`here <https://dl.dropboxusercontent.com/u/58918851/script_info.pdf>`_,
and a discussion about the different implemented ICVs can be found
`here <https://dl.dropboxusercontent.com/u/58918851/icv_info.pdf>`_.

Please, do not hesitate to send a mail to victor.gil.sepulveda@gmail.com
with your questions, criticisms and whatever you think it is not working
or can be done better. Any contribution can help to improve this
software!

TODO
====

-  To improve this documentation (better explanations, more examples and
   downloadable scripts).

-  Refactoring and general improvements:

   -  Total refactoring (Clustering and Clusters are inmutable, hold a
      reference to the matrix -> prototypes are always updated)
   -  Rename script stuff
   -  Rename functions and vars
   -  Minimizing dependencies with scipy
   -  Minimizing dependencies with prody (create my own reader)
   -  Adding its own Hierarchical clustering code (educational
      motivations)
   -  Improve spectral algoritm (add more tests - comparisons with other
      implementations, adding new types)
   -  Improve MPI load balance (i.e. parameter generation must be
      processed in parallel)
   -  Improve test coverage
   -  The script must accept numbers and percentages
   -  Use JSON schema to validate the script. Try to delegate the full
      responsibility of validating to pyProCT (instead of the GUI) [x] -
      Users must be able to comment their scripts (with '//' for
      instance).
   -  When loading a dcd file, we only want to load atomic data of the
      the associated pdb.
   -  Change "compression" by "redundancy\_elimination"
   -  Allow to load all files (or glob) from a folder.

-  Symetry handling:

   -  Symmetry handling for fitting coordinates.
   -  Improve symmetry handling for calculation coordinates (e.g.
      ligands). [x] - Simple chain mapping feature.

-  New algorithms:

   -  Modularity-based (Newman J. 2003)
   -  Passing messages (Frey and Dueck 2007)
   -  Flow simmulation (Stijin van Dongen)
   -  Fuzzy Clustering
   -  `Jarvis-Patrick
      Algorithm <http://www.improvedoutcomes.com/docs/WebSiteDocs/Clustering/Jarvis-Patrick_Clustering_Overview.htm>`_
   -  Others (adaptative spectral clustering flavours)

-  New quality functions.

   -  Balancedness: The sizes of the clusters must be balanced.
   -  J quality function: Cai Xiaoyan Proceedings of the 27th Chinese
      Control Conference
   -  Metastability function (Q) in Chodera et al. J. Chem. Phys. 126
      155101 2007 .
   -  Improve separation quality functions.
   -  New standard separation ICVs (require inmutable prototypes)

          Separation, the clusters themselves should be widely spaced.
          There are three common approaches measuring the distance
          between two different clusters:

          -  Single linkage: It measures the distance between the
             closest members of the clusters.
          -  Complete linkage: It measures the distance between the most
             distant members.
          -  Comparison of centroids: It measures the distance between
             the centers of the clusters.

-  New features:

   -  Refine noise in DBSCAN
   -  Refine a preselected cluster (e.g "noise"), or "heterogeneous".

-  New postprocessing options:

   -  Refinement
   -  Kinetic analysis


