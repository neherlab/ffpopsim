/*
 * main.h
 *
 *  Created on: Dec 28, 2011
 *  Author: Fabio Zanini
 */
#ifndef LOWD_H_
#define LOWD_H_

/* Include directives */
#include <string>
#include "ffpopsim_lowd.h"
#define LOWD_BADARG -1354341
#define NOTHING 1e-10

/* Be verbose? */
#define LOWD_VERBOSE 1

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
int pop_observables();


#endif /* LOWD_H_ */
