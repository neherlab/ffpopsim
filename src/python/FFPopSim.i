/**
* @file FFPopSim.i
* @brief Python 2 bindings of the FFPopSim library.
* @author Richard Neher, Fabio Zanini
* @version 
* @date 2012-04-24
*/
%define DOCSTRING
"C++/Python library for population genetics.

This library offers *two* simulation packages for population genetics: one for
low-dimensional simulations (up to ~15 loci) and one for high-dimensional ones.

Each package is based on a big class that represents a population:

   - ``haploid_lowd`` for low-dimensional populations
   - ``haploid_highd`` for high-dimensional simulations

A simple example routine is the following::

    #####################################
    #   EXAMPLE SCRIPT                  #
    #####################################
    import numpy as np
    import matplotlib.pyplot as plt
    import FFPopSim as h
    
    c = h.haploid_lowd(4)
    c.set_allele_frequencies([0,0.3,0.6,0.9], N=1000) 
    c.evolve(100)
    c.plot_diversity_histogram()
    plt.show()
    #####################################

which evolves a population with 4 loci for 100 generations starting from fixed
allele frequencies, under neutral conditions, and plots the diversity
histogram afterwards.

For more usage examples, please consult the ``tests`` and ``examples`` folders. 
"
%enddef
%module(docstring=DOCSTRING) FFPopSim;

/* Include in the wrap code (note that the multiple header are redundant) */
%{
#define SWIG_FILE_WITH_INIT
#include "../ffpopsim_generic.h"
#include "../ffpopsim_lowd.h"
#include "../ffpopsim_highd.h"
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
/* ffpopsim.h (GENERAL OBJECTS) */
%include "ffpopsim_generic.i";
%include "../ffpopsim_generic.h";

/* ffpopsim_lowd.h (LOW DIMENSIONAL OBJECTS) */
%include "ffpopsim_lowd.i";
%include "../ffpopsim_lowd.h";

/* ffpopsim_highd.h (HIGH DIMENSIONAL OBJECTS) */
%include "ffpopsim_highd.i";
%include "../ffpopsim_highd.h";

/* hivpopulation.h (HIV-SPECIFIC) */
%include "hivpopulation.i";
%include "../hivpopulation.h";
