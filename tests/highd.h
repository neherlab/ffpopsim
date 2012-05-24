/**
 * @file highd.cpp
 * @brief Tests for the high-dimensional simulation library (header).
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-20
 */
#ifndef HIGHD_H_
#define HIGHD_H_

/* Include directives */
#include "hivpopulation.h"
#define HIGHD_BADARG -1354341
#define NOTHING 1e-10

/* Be verbose? */
#define HIGHD_VERBOSE 1

/* main */
int main(int argc, char **argv);

/* generic tests */
int library_access();
int sample_initialize();

/* hypercube testing */
int hc_initialize();
int hc_setting();

/* population testing */
int pop_initialize();
int pop_evolve();
int pop_sampling();
int pop_Hamming();
int pop_divdiv();
int pop_histograms();

/* HIV subclass testing */
int hiv_initialize();
int hiv_evolve();
int hiv_multiple_evolution();

#endif /* HIGHD_H_ */
