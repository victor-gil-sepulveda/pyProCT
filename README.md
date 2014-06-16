pyProCT
==========

pyProCT is an open source cluster analysis software especially adapted for jobs related with structural proteomics. Its approach allows users to define a clustering goal (clustering hypothesis) based on their domain knowledge. This hypothesis will guide the software in order to find the best algorithm and parameters (including the number of clusters) to obtain the result that fulfills their expectatives. In this way users do not need to use cluster analysis algorithms as a black box, improving this way the results.
pyProCT not only generates a resulting clustering, it also implements some use cases like the extraction of representatives or trajectory redundance elimination.

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
		 "base": /home/john/ClusteringProject"
	}
}
```
This is an example of "global" 

### Data

### Clustering

### Postprocessing

### Checking the script
As the "script" is indeed a JSON object, any JSON validator can be used to discover the errors in case of script loading problems. A good example of such validators is [JSONLint](http://jsonlint.com/). 

### Learning more
A more detailed explanation of the script contents can be found [here](pdf/script_info.pdf), and a discussion about the different implemented ICVs can be found [here](pdf/icv_info.pdf). 

## Using pyProCT as part of other programs 
pyProCT has been written as a collection of modules coupled by means of the different handlers and 
Algorithms can be used separately if the correct for example:


The necessary documentation to use pyProCT classes is written inside the code. It has been extracted [here]() and [here](). We are currently trying to improve this documentation with better explanations and examples. 

### Using it as a separate program from other Python script



## A parallel execution example


# Documentation 

We are still experimenting to see which documentation generator fits better with us. Currently we have two versions of the documentations: one using [Sphinx](http://sphinx-doc.org/) and the other using [Doxygen](http://www.stack.nl/~dimitri/doxygen/)+[doxpy](http://code.foosel.org/doxypy). See them [here](pyproct/docs/_build/html/index.html) and [here](pyproct/docs/doxyxml/html/index.html).










