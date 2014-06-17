pyProCT
==========

pyProCT is an open source cluster analysis software especially adapted for jobs related with structural proteomics. Its approach allows users to define a clustering goal (clustering hypothesis) based on their domain knowledge. This hypothesis will guide the software in order to find the best algorithm and parameters (including the number of clusters) to obtain the result that fulfills their expectatives. In this way users do not need to use cluster analysis algorithms as a black box, improving this way the results.
pyProCT not only generates a resulting clustering, it also implements some use cases like the extraction of representatives or trajectory redundance elimination.

[Table of Contents**](http://doctoc.herokuapp.com/)

- [pyProCT](#user-content-pyproct)
	- [Installation](#user-content-installation)
	- [Using pyProCT as a standalone program](#user-content-using-pyproct-as-a-standalone-program)
		- [Global](#user-content-global)
		- [Data](#user-content-data)
		- [Clustering](#user-content-clustering)
			- [generation](#user-content-generation)
			- [algorithms](#user-content-algorithms)
			- [evaluation](#user-content-evaluation)
		- [Postprocessing](#user-content-postprocessing)
		- [Checking the script](#user-content-checking-the-script)
		- [Learn more](#user-content-learn-more)
	- [Using pyProCT as part of other programs](#user-content-using-pyproct-as-part-of-other-programs)
		- [Using it as a separate program from other Python script](#user-content-using-it-as-a-separate-program-from-other-python-script)
	- [A parallel execution example](#user-content-a-parallel-execution-example)
- [Documentation](#user-content-documentation)

## Installation
pyProCT is quite easy to install using *pip*. Just write:

```Shell
> sudo pip install pyProCT
```
And *pip* will take care of all the dependencies (shown below).

<img src="img/dependencies.png"> </img>

- <img src="img/warning.png"></img>  It is recommended to install Numpy and Scipy before starting the installation using your OS software manager. You can try to download and install them [manually](http://docs.scipy.org/doc/numpy/user/install.html) if you dare.

- <img src="img/warning.png"></img>  mpi4py is pyProCT's last dependency. It can give problems when installing it in OS such as SUSE. If the installation of this last package is not succesful, pyProCT can still work in Serial and Parallel (using *multiprocessing*) modes.

## Using pyProCT as a standalone program

The preferred way to use pyProCT is throug a JSON "script" that describes the clustering task. It can be executed using the following line in your shell:

```Shell
> python -m pyproct.main script.json
```

The JSON script has 4 main parts, each one dealing with a different aspect of the clustering pipeline. This sections are:
* global: Handles workspace and scheduler parameterization.
* data: Handles distance matrix parameterization.
* clustering: Handles algorithms and evaluation parameterization.
* preprocessing: Handles what to do with the clustering we have calculated.

```JSON
{
	"global":{},
	"data":{},
	"clustering":{},
	"postprocessing":{}
}
```

### Global
```JSON
{
	"control": {
		"scheduler_type": "Process/Parallel",
		"number_of_processes": 4
	},
	"workspace": {
		 "tmp": "tmp",
		 "matrix": "matrix",
		 "clusterings": "clusterings",
		 "results": "results",
		 "base": "/home/john/ClusteringProject"
	}
}
```
This is an example of _global_ section. It describes the work environment (workspace) and the type of scheduler that will be built. Defining the subfolders of the wokspace is not mandatory, however it may be convenient in some scenarios (for instance, in serial multiple clustering projects, sharing the _tmp_ folder would lower the disk usage as at each step it will be ooverwritten).

This is a valid global section using a serial scheduler and default names for workspace inner folders:
```JSON
{
	"control": {
		"scheduler_type": "Serial"
	},
	"workspace": {
		 "base": "/home/john/ClusteringProject"
	}
}
```
pyProCT allows the use of 3 different schedulers that help to improve the overall performance of the software by parallelizing some parts of the code. The available schedulers are "Serial", "Process/Parallel" (uses Python's [multiprocessing](https://docs.python.org/2/library/multiprocessing.html)) and "MPI/Parallel" (uses MPI through the module [mpi4py](http://mpi4py.scipy.org/)).

### Data
The _data_ section defines how pyProCT must build the distance matrix that will be used by the compression algorithms. Currently pyProCT offers up to three options to build that matrix: "load", "rmsd" and "distance"
- rmsd: Calculates a all vs all rmsd matrix using any of the [pyRMSD](https://github.com/victor-gil-sepulveda/pyRMSD) calculators available. It can calculate the RMSD of the fitted region (defined by [Prody](http://prody.csb.pitt.edu/) compatible selection string in _fit_selection_) or one can use one selection to superimpose and another to calculate the rmsd (_calc_selection_) .
-  distance: After superimposing the selected region it calculates the all vs all distances of the geometrical center of the region of interest (_body_selection_).
- load: Loads a precalculated matrix.

```JSON
{
	"type": "pdb_ensemble",
	"files": [
		"A.pdb",
		"B.pdb"
	],
	"matrix": {
		"method": "rmsd",
		"parameters": {
			"calculator_type": "QCP_OMP_CALCULATOR",
			"fit_selection": "backbone",
		}
		"image": {
			"filename": "matrix_plot"
		},
		"filename":"matrix"
	}
}
```
The matrix can be stored if the _filename_ property is defined. The matrix can also be stored as an image if the _image_ property is defined.

pyProCT can currently load _pdb_ and _dcd_ files. When using _pdb_ files, files can be loaded in two ways:
1. Using a list of file paths.
2. Using a list of file objects:
```JSON
{ 
	"file": ... , 
	"base_selection": ... 
}
```
Where _base_selection_ is a [Prody](http://prody.csb.pitt.edu/) compatible selection string. Loading files this way can help in cases where not all files have structure with the same number of atoms: _base_selection_ should define the common region between them (if a 1 to 1 map does not exist, the RMSD calculation will be wrong).

3. Only for _dcd_ files:
```JSON
{ 
	"file": ...,
	"atoms_file": ..., 
	"base_selection": ... 
}
```
Where _atoms_file_ is a _pdb_ file with at least one frame that holds the atomic information needed by the _dcd_ file.

### Clustering
The _clustering_ section is divided in 3 other subsections: 

```JSON
{
	"generation": {
		“method”: “generate”
	},
	"algorithms": {
		...
	},
	"evaluation": {
		...
	}
}
```

#### generation 
Defines how the clustering will be generated (_load_ or _generate_). if _load_ is chosen, the section must contain the clustering that may be used. Ex:

```JSON
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
```

#### algorithms
If pyProCT has to generate the clustering, this section defines the algorithms that will be used as well as their parameters (if necessary). The currently available algorithms are : _kmedoids_, _hierarchical_, _dbscan_, _gromos_, _spectral_ and _random_. Each algorithm can store its list of parameters, however the preferred way to work with pyProCT is to let it automatically generate them. Almost all algorithms accept the property _max_, that defines the maximum amount of parameter collections that will be generated for that algorithm.

Ex.

```JSON
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
```

Algorithm parameters can be explicitly written:

```JSON
"kmedoids": {
	"seeding_type": "RANDOM",
	"max": 50,
	"tries": 5,
	"parameters":[{"k":4},{"k":5},{"k":6}]
}
```

#### evaluation
This section holds the _Clustering Hypothesis_, the core of pyProCT. Here the user can define how will the expected cluste

Ex.
```JSON
{
	"maximum_noise": 15,
	"minimum_cluster_size": 50,
	"maximum_clusters": 200,
	"minimum_clusters": 6,
	"evaluation_criteria": { 
		...
	},          
	"query_types": [ ... ]          
}
```
### Postprocessing

### Checking the script
As the "script" is indeed a JSON object, any JSON validator can be used to discover the errors in case of script loading problems. A good example of such validators is [JSONLint](http://jsonlint.com/). 

### Learn more
A more detailed explanation of the script contents can be found [here](pdf/script_info.pdf), and a discussion about the different implemented ICVs can be found [here](pdf/icv_info.pdf).

## Using pyProCT as part of other programs 
pyProCT has been written as a collection of modules coupled by means of the different handlers and 
Algorithms can be used separately if the correct for example:

The necessary documentation to use pyProCT classes is written inside the code. It has been extracted [here]() and [here](). We are currently trying to improve this documentation with better explanations and examples. 

### Using it as a separate program from other Python script

Loading results
generating scripts programatically

## A parallel execution example

# Documentation 
We are still experimenting to see which documentation generator fits better with us. Currently we have two versions of the documentations: one using [Sphinx](http://sphinx-doc.org/) and the other using [Doxygen](http://www.stack.nl/~dimitri/doxygen/)+[doxpy](http://code.foosel.org/doxypy). See them [here](pyproct/docs/_build/html/index.html) and [here](pyproct/docs/doxyxml/html/index.html). We will possibly publish it in a cloud solution like [readthedocs.org](https://readthedocs.org/)
