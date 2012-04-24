/*
author:     Fabio Zanini
date:       17 September 2011
content:    Extend Python with a C++ function using Richard's evolution library.
*/
#ifndef _EVOLUTION_H
#define _EVOLUTION_H


int
check_normalization(int gts, double* frequencies);

/*
N.B.: the last input arguments are modified in-place via SWIG magic.
      They are given in Python via a 1D, 2**L-long float numpy array.
*/
int 
evolve(double N, int generations, int gts, int loci, double* hoprates, int gts2, double* frequencies);
//	    ^		     ^	    ^		^		^	  ^	    	  ^
//	    |		     |	    |___________|_______________|	  |_______________|
//	    |		     |		        |				  |
// population size  number of generations       |				  |
//  					        |				  |
//			2^L x L (flattened) matrix for hopping rates   2^L array for the distribution
//					SWIG INPUT				SWIG INPLACE


int
evolve_consensus(double N, int generations, double mutrate, int loci, double* fitness_additive, int gts2, double* frequencies, int reps, long seed=0);
//			^	     ^	  	  ^		^	   		^	   ^	    	     ^		  ^		^
//			|	     |	  	  |		|_______________________|	   |_________________|		  |		|
//			|	     |		  |		   	|				  |		  number of repetitions	|
//	 population size  number of generations   |		   	|				  |		  for each time point	|
//	  					  |		   	|				  |					|
//						  |     additive fitness coefficients         2^L array for the distribution		seed for random generator
//						mutation	SWIG INPUT				SWIG INPLACE
//						rate


#endif // _EVOLUTION_H
