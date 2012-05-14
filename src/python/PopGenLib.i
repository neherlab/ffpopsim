/**
* @file hivpopulation.i
* @brief Python 2 bindings of the hivpopulation class.
* @author Richard Neher, Boris Shraiman, Fabio Zanini
* @version 
* @date 2012-04-24
*/
%module(docstring="PopGenLib library for population genetics.") PopGenLib;
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

/**************************************************************
 * CODE TO BE WRAPPED
 *************************************************************/

/* popgen.h (GENERAL OBJECTS) */
%include "popgen.i"
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
