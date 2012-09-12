/**
 * @file lowd.cpp
 * @brief Tests for the low-dimensional simulation library.
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-20
 */
/* Include directives */
#include <string>
#include "ffpopsim_lowd.h"
#define LOWD_BADARG -1354341
#define NOTHING 1e-10

/* Be verbose? */
#define LOWD_VERBOSE 1

/* Test generic library access */
int library_access() {
	index_value_pair_t ivp;
	ivp.index = 4;
	ivp.val = 4.5;

	if(LOWD_VERBOSE) {
		cerr<<"index_value_pair: index = "<<ivp.index;
		cerr<<", val = "<<ivp.val<<endl;
	}
	return 0;
}


/* Test sample initialization */
int sample_initialize() {
	int N = 5;

	sample sam;
	sam.set_up(N);

	for(int i=0; i< N; i++)
		sam.values[i] = i + 3;
	if(LOWD_VERBOSE)
		cerr<<"Mean of sample: "<<sam.calc_mean()<<endl;
	return 0;
}

/* Test hypercube_lowd initialization */
int hc_initialize() {
	hypercube_lowd hc;
	int status = hc.set_up(4,3);
	if(LOWD_VERBOSE){
		cerr<<"Memory allocation = "<<status<<endl;
		cerr<<"Dimension: "<<hc.get_dim()<<endl;
	}
	return 0;
}

/* Test setting coefficients on hypercube_lowd */
int hc_setting() {
	int L = 2;
	double additive[L];

	for(int i=0; i<L;i++) {
		additive[i] = 3;
	}

	hypercube_lowd hc(3, L);
	hc.reset();
	hc.additive(additive);
	hc.fft_coeff_to_func();
	
	if(LOWD_VERBOSE){
		cerr<<"Func values: ";
		for(int gt=0; gt < (1<<L); gt++)
			cerr<<hc.get_func(gt)<<" ";
		cerr<<endl;

		cerr<<"Coeff values: ";
		for(int gt=0; gt < (1<<L); gt++)
			cerr<<hc.get_coeff(gt)<<" ";
		cerr<<endl;
	}
	return 0;
}

/* Test population initialization */
int pop_initialize() {
	int L = 4;

	haploid_lowd pop(L, 3);
	if(LOWD_VERBOSE)
		cerr<<"L = "<<pop.L()<<endl;	
	return 0;	
}

/* Test evolution from allele frequencies */
int pop_evolve_af() {
	int L = 4;
	int N = 10000;

	haploid_lowd pop(L);

	double freq[1<<L];
	for(int i=0; i<(1<<L);i++)
		freq[i] = 1.0 / (1<<L);

	// Set allele frequencies
	pop.set_allele_frequencies(freq, N);

	double* rr = new double[L-1];
	for(int i=0; i < L-1; i++)
		rr[i] = 0.01;
	pop.set_recombination_rates(rr);	
	pop.set_mutation_rates(1e-2);
	pop.evolve(5);

	if(LOWD_VERBOSE) {
		cerr<<"Population size: "<<pop.N()<<endl;
	}
	return 0;	
}

/* Test evolution from genotype frequencies */
int pop_evolve_gf() {
	int L = 4;
	int N = 10000;

	haploid_lowd pop(L);

	// Start with wildtype and a single mutant
	index_value_pair_t ivp(0, N/2);
	vector<index_value_pair_t> gts;
	gts.push_back(ivp);
	ivp.index = 1;
	gts.push_back(ivp);
	pop.set_genotypes(gts);

	double* rr = new double[L-1];
	for(int i=0; i < L-1; i++)
		rr[i] = 0.01;
	pop.set_recombination_rates(rr);	
	pop.set_mutation_rates(1e-2);
	pop.evolve(5);

	if(LOWD_VERBOSE) {
		cerr<<"Population size: "<<pop.N()<<endl;
	}
	return 0;	
}


/* Test genotype and allele frequency */
int pop_observables() {
	int L = 4;
	int N = 100;

	haploid_lowd pop(L,3);
	double freq[1<<L];
	for(int i=0; i<(1<<L);i++)
		freq[i] = 1.0 / (1<<L);
	pop.set_allele_frequencies(freq, N);
	pop.set_mutation_rates(1e-2);
	pop.evolve(5);

	if(LOWD_VERBOSE) {
		cerr<<"Allele freq[0]: "<<pop.get_allele_frequency(0)<<endl;
		cerr<<"Fitness mean and variance: "<<pop.get_fitness_statistics().mean<<", "<<pop.get_fitness_statistics().variance<<endl;

	}

	return 0;	
}

/* MAIN */
int main(int argc, char **argv) {
	int status;
	if (argc > 1) {
		cout<<"Usage: "<<argv[0]<<endl;
		status = 1;
	} else {
		status = library_access();
		status += sample_initialize();
		status += hc_initialize();
		status += hc_setting();
		status += pop_initialize();
		status += pop_evolve_af();
		status += pop_evolve_gf();
		status += pop_observables();
	}
	cout<<"Number of errors: "<<status<<endl;
	return status;
}
