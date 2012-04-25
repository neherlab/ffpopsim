/**
* @file hivpopulation.i
* @brief Python 2 bindings of the hivpopulation class.
* @author Richard Neher, Boris Shraiman, Fabio Zanini
* @version 
* @date 2012-04-24
*/
%module(docstring="Model for an HIV population.") hivpython
/* Include in the wrap code */
%{
#define SWIG_FILE_WITH_INIT
#include "hivpython.h"
%}

%ignore check_input;

/* Numpy magic to output arrays */
%include "numpy.i"
%init %{
import_array();
%}

/*%apply (int DIM1, double* INPLACE_ARRAY1) {(int gts2, double* frequencies)};
%apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int gts, int loci, double* hoprates)};
%apply (int DIM1, double* IN_ARRAY1) {(int loci, double* fitness_additive)};
*/

%include "hivpython.h"

