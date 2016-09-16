---
title: 'ecopath_matlab: A Matlab-based implementation of the Ecopath food web algorithm'
tags:
  - Ecopath
  - food web
  - network analysis
authors:
 - name: Kelly A. Kearney
   affiliation: University of Washington, JISAO
date: 14 September 2016
bibliography: paper.bib
---

# Summary

The Ecopath with Ecosim model [@christensen2004] is one of the most widely-used ecosystem modeling frameworks in fisheries science, and its use has contributed to over 400 scientific publications [@Akoglu2015].  The software is centered on a mass-balance algorithm (Ecopath) that uses commonly-measured metrics of fish and fisheries to quantify key features of a food web network, including biomass/energy of all state variables in an ecosystem (living groups, detrital pools, and fishing fleets) and the fluxes between them.

While the Ecopath with Ecosim software is freely available at www.ecopath.org, and its source code is available by request from the authors for development of external plugins [@Steenbeek2015], several design factors of that software limit its ease of use for extended development.  First, it is currently written for the Microsoft .NET framework, a language that is used by very few researchers in the fields of fisheries and ocean sciences, and is difficult to use outside of Microsoft operating systems.  Second, the code for the primary scientific algorithms is intertwined with the code for the graphical user interface within which those algorithms are intended to be used.  This limits the ease with which the Ecopath algorithm can be used outside of the graphical user interface, including automation of inputs, extended analysis of outputs, and coupling of the algorithm to other models.  

The desire to accomplish extended Ecopath-based research, such as ensemble model generation or coupling of food web models to ocean biogeochemical models, has led to several concurrent and independent efforts to translate the core Ecopath algorithm into languages more familiar to the fisheries and oceanographic communities.  This includes Rpath [@Lucey] for R, a language favored for statistical analysis in fisheries science, and EwE-F [@Akoglu2015] for Fortran, targeted towards use with hydrodynamic and biogeochemical models.

This package adds a Matlab-based version of Ecopath to this collection.  The Matlab language is commonly used in ocean and atmospheric sciences, primarily for data analysis and vizualization, as well as for preliminary model developement and testing.  The goal of this implementation is to provide a lightweight replication of Ecopath, allowing for easy automation of input and output and freeing that input and output to be used with the wide variety of existing Matlab functions, toolboxes, and graphical capabilities.  The primary features of this implementation are:

- Full replication of the Ecopath algorithm, including support for multi-stanza groups and support for less common input options (e.g. groups with diet import, partial primary producers, partial use of habitat area, etc.).
- Routines to explore model parameter uncertainty through the use of model ensembles.
- Routines to calculate common network analysis metrics, as well as to convert an Ecopath model to a Matlab graph object for further graph theory-style analysis.
- Functions to import and export data between the existing flavors of Ecopath, including Ecopath with Ecosim database files, Rpath comma-delimited files, and EwE-F tab-delimited text files.

# Publications

This code, or previous versions of it, has been used in the following publications:

- [@Kearney2012]
- [@Kearney2012a]
- [@Kearney2013]
- [@Kearney2015]
- [@Guesnet2015]

# References
