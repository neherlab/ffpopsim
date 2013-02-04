/**
 * Copyright (c) 2012-2013, Richard Neher, Fabio Zanini
 * All rights reserved.
 *
 * This file is part of FFPopSim.
 *
 * FFPopSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FFPopSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FFPopSim. If not, see <http://www.gnu.org/licenses/>.
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

/* STL typemaps */
%include <typemaps.i>
%include <std_string.i>
%include <std_vector.i>
%include <std_list.i>
%include <std_map.i>

/* Numpy magic to output arrays */
/*%include "pyfragments.swg"*/
%include "numpy.i";
%init %{
import_array();
%}

/* STL typemaps */
%include <typemaps.i>
%include <std_vector.i>
%template(_intVector) std::vector<int>;
/*%template(_doubleVector) std::vector<double>;*/
%template(vector_tree_step) std::vector<step_t>;
%template(vector_tree_key) std::vector<tree_key_t>;
%template(list_tree_key) std::list<tree_key_t>;
%template(map_key_edge) std::map<tree_key_t, edge_t>;
%template(map_key_node) std::map<tree_key_t, node_t>;
%template(vector_polymorphism) std::vector<poly_t>;

/* FFPopSim typemaps */
%include "ffpopsim_typemaps.i";

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
