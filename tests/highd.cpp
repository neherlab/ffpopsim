/**
 * @file highd.cpp
 * @brief Tests for the high-dimensional simulation library.
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-20
 */
#include "highd.h"

/* MAIN */
int main(int argc, char **argv){
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
		status += pop_evolve();
		status += pop_observables();
	}
	cout<<"Number of errors: "<<status<<endl;
	return status;
}


/* Test generic library access */
int library_access() {
	index_value_pair ivp;
	ivp.index = 4;
	ivp.val = 4.5;

	if(HIGHD_VERBOSE) {
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
	if(HIGHD_VERBOSE)
		cerr<<"Mean of sample: "<<sam.calc_mean()<<endl;
	return 0;
}

/* Test hypercube initialization */
int hc_initialize() {
	hypercube_function hc;
	int status = hc.set_up(4,3);
	if(HIGHD_VERBOSE){
		cerr<<"First constructor: Memory allocation = "<<status<<endl;
		cerr<<"First constructor: Dimension: "<<hc.get_dim()<<endl;
	}
	hypercube_function hc2(4,3);
	if(HIGHD_VERBOSE){
		cerr<<"Second constructor: Dimension: "<<hc2.get_dim()<<endl;
	}
	return 0;
}

/* Test setting coefficients on hypercube */
int hc_setting() {
	int L = 1000;

	hypercube_function hc(L,3);

	double values[] = {0.1, 0.043, 1.3, -4.3};
	int myints[] = {16,2,77,29};
	vector<int> loci;
	for(int i=0; i<4; i++) {
		loci.assign(1, myints[i]);
		hc.add_coefficient(values[i], loci);
		loci.clear();
	}	
	
	if(HIGHD_VERBOSE){
		cerr<<"Coefficient indices and values: ";
		for(int i=0; i<4; i++) {
			cerr<<"("<<myints[i]<<","<<hc.get_additive_coefficient(myints[i])<<") ";
		}
		cerr<<endl;
	}
	return 0;
}

/* Test population initialization */
int pop_initialize() {
	int L = 1000;
	int N = 100;

	haploid_clone pop(N, L, 3, 1);
	if(HIGHD_VERBOSE)
		cerr<<"L = "<<pop.L()<<", N = "<<pop.get_target_pop_size()<<endl;	
	return 0;	
}

/* Test evolution */
int pop_evolve() {
	int L = 1000;
	int N = 100;

	haploid_clone pop(N, L, 3, 1);

	pop.set_mutation_rate(1e-2);
	pop.set_crossover_rate(1e-2);
	pop.init_genotypes();		// start with a population of the right size
	pop.evolve(5);

	return 0;	
}

/* Test genotype and allele frequency */
int pop_observables() {
	int L = 1000;
	int N = 100;

	haploid_clone pop(N, L);

	pop.set_mutation_rate(1e-3);
	pop.set_outcrossing_probability(1e-2);
	pop.set_crossover_rate(1e-2);
	pop.set_recombination_model(CROSSOVERS);
	pop.init_genotypes();		// start with a population of the right size

	vector <int> loci;
	for(int i=0; i< 100; i++) {
		loci.assign(1, i);
		pop.add_fitness_coefficient(0.1, loci);
		loci.clear();
	}

	pop.calc_everybodies_traits();
	pop.evolve(10);
	pop.calc_stat();


	if(HIGHD_VERBOSE) {
		double af = 0;
		double tmp;
		for(int i=0; i< L; i++) {
			tmp = pop.get_allele_frequency(i);
			af = MAX(af, tmp);
		}

		cerr<<"Generation: "<<pop.get_generation()<<endl;
		cerr<<"Max allele freq: "<<af<<endl;
		cerr<<"Fitness mean and variance: "<<pop.fitness_mean()<<", "<<pop.fitness_var()<<endl;
	}

	return 0;	
}


