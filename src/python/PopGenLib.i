/**
* @file hivpopulation.i
* @brief Python 2 bindings of the hivpopulation class.
* @author Richard Neher, Boris Shraiman, Fabio Zanini
* @version 
* @date 2012-04-24
*/
%module(docstring="PopGenLib library for population genetics.") PopGenLib
/* Include in the wrap code (note that the multiple header are redundant) */
%{
#define SWIG_FILE_WITH_INIT
#include "../popgen.h"
#include "../popgen_highd.h"
#include "../hivpopulation.h"
%}

/* kwargs support (does not work well with overloaded functions) */
%feature("kwargs");

/* Numpy magic to output arrays */
/*%include "pyfragments.swg"*/
%include "numpy.i"
%init %{
import_array();
%}


/**************************************************************
 * CODE TO BE WRAPPED
 *************************************************************/

/* popgen.h (GENERAL OBJECTS) */
%rename(index_value_pair_t) index_value_pair;
%rename(stat) stat_t;
%ignore SAMPLE_ERROR;
%ignore sample;
%include "../popgen.h"

/* popgen_highd.h (HIGH DIMENSIONAL OBJECTS) */
%include "popgen_highd.i"
%include "../popgen_highd.h"

/* hivpopulation.h (HIV-SPECIFIC) */
%include "hivpopulation.i"
%include "../hivpopulation.h"
