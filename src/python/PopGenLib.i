/**
* @file hivpopulation.i
* @brief Python 2 bindings of the hivpopulation class.
* @author Richard Neher, Boris Shraiman, Fabio Zanini
* @version 
* @date 2012-04-24
*/
%define DOCSTRING
"PopGenLib library for population genetics.

This library offers two simulation packages for population genetics. Each is controlled by a basic class:

- haploid_lowd: low-dimensional populations (genomes shorter than ~20 loci)
- haploid_clone:  high-dimensional simulations (genomes longer than ~20 loci)

The library is written in C++ and offers a Python interface. A simple evolution routine could be the following:

#####################################
#   EXAMPLE SCRIPT                  #
#####################################
import numpy as np
import matplotlib.pyplot as plt
import PopGenLib as h

c = h.haploid_lowd(4, 1000)
c.init_frequencies([0,0.3,0.6,0.9]) 
c.evolve(10)
c.plot_fitness_histogram()
plt.show()
#####################################

which evolves a 4-loci population for 10 generations starting from fixed allele frequencies, and plots the
fitness histogram afterwards. Please look into the 'tests' folder for more usage examples. 

Requirements:
- numerical Python is used extensively in this library, and is *strongly* recommended for all users. We suggest to import numpy explicitely before using the library (but it will work in any case).
- matplotlib is used in the plot functions. As long as you do not call those, you can live without it. However, we suggest to import it explicitely before using the library.

*Note*: the Python interface does not offer the full functionality of the underlying library, nor as much
flexibility. If you need to perform a peculiar type of simulations that is not accessible by this Python
interface, consider C++ (or Python) subclassing.
"
%enddef
%module(docstring=DOCSTRING) PopGenLib;

/* Include in the wrap code (note that the multiple header are redundant) */
%{
#define SWIG_FILE_WITH_INIT
#include "../popgen.h"
#include "../popgen_lowd.h"
#include "../popgen_highd.h"
#include "../hivpopulation.h"
%}

/* Numpy magic to output arrays */
/*%include "pyfragments.swg"*/
%include "numpy.i";
%init %{
import_array();
%}

/* activate autodoc */
%feature("autodoc", "1");

/**************************************************************
 * CODE TO BE WRAPPED
 *************************************************************/
/* popgen.h (GENERAL OBJECTS) */
%include "popgen.i";
%include "../popgen.h";

/* popgen_lowd.h (LOW DIMENSIONAL OBJECTS) */
%include "popgen_lowd.i";
%include "../popgen_lowd.h";

/* popgen_highd.h (HIGH DIMENSIONAL OBJECTS) */
%include "popgen_highd.i";
%include "../popgen_highd.h";

/* hivpopulation.h (HIV-SPECIFIC) */
%include "hivpopulation.i";
%include "../hivpopulation.h";
