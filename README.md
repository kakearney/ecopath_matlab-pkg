<img src="./matlabEcopathLogo.png" width="100">

# A Matlab implementation of Ecopath 

This package provides a Matlab-based implementation of the Ecopath algorithm (part of the popular fisheries ecosystem modeling tool [Ecopath with Ecosim](http://www.ecopath.org)), as well as a few functions for further analysis and visualization of an Ecopath-style food web.

This software is intended for people already familiar with the Ecopath concept who wish to perform that particular calculation, as well as related analysis of food web properties, outside of the GUI environment provided by the original EwE software.  It assumes a basic working knowledge of Matlab.

## Getting Started

### Prerequisites

This software requires [Matlab R2015b](http://www.mathworks.com/products/matlab/) or later.

For full use, it also requires the [Statistics and Machine Learning Toolbox](http://www.mathworks.com/products/statistics/).

No OS limitations beyond those required for Matlab itself.
  
### Downloading

Git users can clone directly:

`git clone git@github.com:kakearney/ecopath_matlab-pkg`

Alternatively, you may download a zipped version of the source code via the _Clone or Download_ > _Download zip_ button above, or from the [ecopath_matlab](http://www.mathworks.com/matlabcentral/fileexchange/40082-ecopathlite-m--a-matlab-implementation-of-ecopath "FEX ecopath_matlab") entry on the MatlabCentral File Exchange.  The File Exchange entry is updated daily from this repository.

### Installation

Once downloaded (and unzipped, if necessary), all code subfolders in the main folder (i.e all folders except the `examples` and `html` folders) need to be added to your Matlab path:

```matlab
% Replace with location of the downloaded folder
pth = './ecopath_matlab-pkg/';

addpath(fullfile(pth, 'ConsoleProgressBar'); 
addpath(fullfile(pth, 'aggregate'); 
addpath(fullfile(pth, 'cellstr2'); 
addpath(fullfile(pth, 'cprintf'); 
addpath(fullfile(pth, 'ecopath_matlab'); 
addpath(fullfile(pth, 'legendflex'); 
addpath(fullfile(pth, 'readtext'); 
addpath(fullfile(pth, 'regexpfound'); 
addpath(fullfile(pth, 'setgetpos_V1.2'); 
addpath(fullfile(pth, 'wraptext'); 
```

## Package Contents

### The `ecopathmodel` class

This package centers around a custom Matlab class, `ecopathmodel`.  The properties of an `ecopathmodel` object hold the typical input parameters associated with a single ecosystem food web; the methods provide functions to calculate Ecopath mass balance.

Constructor Summary

*  `ecopathmodel`:	Create an ecopathmodel object 

Property Summary

*  `dc`:	Table of diet composition data 
*  `df`:	Table of detritus fates 
*  `discard`:	Table of fisheries discards 
*  `discardFate`:	Table of discard fates 
*  `fleet`:	Names corresponding to each fishing gear/fleet in the model 
*  `groupdata`:	Table of group-related parameters 
*  `landing`:	Table of fisheries landings 
*  `name`:	Names corresponding to each group in the model 
*  `ngear`:	Number of fishing gears/fleets in model 
*  `ngroup`:	Number of groups (living and detrital) in model 
*  `nlive`:	Number of living groups (non-detrital) in model 
*  `pedigree`:	Table of pedigree values applied to parameters 
*  `stanza`:	Names corresponding to each multi-stanza set 
*  `stanzadata`:	Table of multi-stanza-set-related parameters 

Method Summary

*  `addpedigree`:	Add entries to the pedigree table 
*  `calcstanza`:	Calculate B and Q/B values for multi-stanza Ecopath groups 
*  `checkstanza`:	Fill in (or validate) B, QB, and BA values for stanzas 
*  `combinegroups`:	Combine groups and/or fleets in an ecopathmodel object 
*  `createensemble`:	Build an ensemble of Ecopath model parameters 
*  `displaybasic`:	Prints ecopath results to screen 
*  `ecopath`:	Rewrite of Ecopath algorithms 
*  `getpedigreevals`:	Extract values corresponding to pedigree entries 
*  `graph`:	Convert ecopathmodel object to a diagraph object 
*  `networkindices`:	Calculate ecological network indices 
*  `sort`:	Sort groups and fleets in an ecopathmodel object 
*  `sortbytrophic`:	Sort ecopathmodel object groups by trophic level 
*  `stanzaindices`:	Extract indices of stanza groups, in order of age 
*  `subpedigreevalues`:	Replace values in ecopathmodel based on pedigree 
*  `unitconvert`:	Convert units of parameter values 

### Additional functions

A few additional functions are provided alonside the `ecopathmodel` class, including functions for import and export of data, as well as a few helper functions that are called by the `ecopathmodel` class methods but can also be called independently by users:

* `ecopathmodel2rpath`: Print ecopathmodel data to comma-delimited files
* `editstanzacalcs`: Replicate multi-stanza calculations from Ecopath
* `mdb2ecopathmodel`: Create ecopathmodel object from EwE6 data file
* `networkindices`: Calculate network indices for a food web
* `rpath2ecopathmodel`: Create ecopathmodel object from Rpath data files
* `trophiclevel`: Estimates trophic level of food web members

## Usage

Documentation for each function and method in this package is provided via standard Matlab function headers, accessed in Matlab via the `help` command.  The reference page for the `ecopathmodel` class (including links to descriptions of all properties and methods) can be accessed by typing the following in the Matlab Command Window:

```matlab
doc ecopathmodel
```

For an overview of the class, along with several examples of usage, please see the overview document: `ecopathmodel_overview.m`.  A published version of this file can be found in [html/ecopathmodel_overview.html](https://rawgit.com/kakearney/ecopath_matlab-pkg/master/html/ecopathmodel_overview.html).




