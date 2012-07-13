/**
 * @file popgen_lowd.h
 * @brief Header file for low-dimensional simulations
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-04-19
 *
 * Copyright (c) 2012, Richard Neher,Fabio Zanini
 * All rights reserved.

 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef FFPOPSIM_LOWD_H_
#define FFPOPSIM_LOWD_H_
#include "ffpopsim_generic.h"

#define HC_MEMERR -131545		//memory error code
#define HC_BADARG -131546		//bad argument error code
#define HC_VERBOSE 0			//debugging: if set to one, each function prints out a message into the error stream
#define HC_FUNC 1			//hypercube_lowd.func is up-to-date
#define HC_COEFF -1			//hypercube_lowd.coeff is up-to-date
#define HC_FUNC_EQ_COEFF 0		//hypercube_lowd.func equal hypercube_lowd.coeff

using namespace std;

/**
 * @brief Binary hypercube_lowd used in low-dimensional simulations.
 *
 * This class is a generic object that can be used to represent various things, e.g.
 * - the fitness landscape or any other phenotypic landscape;
 * - the genotype frequencies of a population with a genome of size L.
 *
 * If you are planning to model a whole population evolving on the hypercube_lowd, see the class haploid_lowd.
 *
 * Notes on scalability:
 * - The number of genotypes to store increases as \f$2^L\f$, where L is the number of sites
 * - The number of recombination intermediates to be stored increases as \f$3^L\f$, this class can thus only be used for \f$L\f$ up to 20 or so.
 * - The population size N is actually unimportant, as far as it can be stored as a long integer. In other words, this class scales with N like O(1).
 */
class hypercube_lowd
{
public:
	//dimension of the hypercube_lowd
	int dim;

	int state;				//takes values HC_FUNC, HC_COEFF, HC_HC_FUNC_EQ_COEFF, depending on the current state of hypercube_lowd
	double *coeff;				//array holding 2^N coefficients: a entry 0101001101 corresponds to a term with spins at each 1
	double *func;				//array holding the values of the function on the hypercube_lowd

	int *order;				//Auxiliary array holding the number of spins, i.e. the number of ones of coeff[k]

	// construction / destruction
	hypercube_lowd();
	hypercube_lowd(int dim_in, int s=0);
	~hypercube_lowd();
	int set_up(int dim_in, int s=0);

	// set coefficients
	int gaussian_coefficients(double* vark, bool add=false);
	int additive(double* additive_effects, bool add=false);
	int init_rand_gauss(double sigma, bool add=false);
	int init_list(vector<index_value_pair_t> iv, bool add=false);
	int init_coeff_list(vector <index_value_pair_t> iv, bool add=false);
	void calc_order();
	void set_state(int s){state=s;}

	//in and out
	int read_coeff(istream &in);
	int write_func(ostream &out);
	int write_coeff(ostream &in,  bool label=false);
	int read_func(istream &out);
	int read_func_labeled(istream &in);

	//analysis
	int signature(int point);

	//transform from coefficients to function and vice versa
	int fft_func_to_coeff();
	int fft_coeff_to_func();

	//read out
	int get_state() {return state;}
	unsigned int get_dim(){return dim;}
	unsigned int get_seed() {return seed;}
	double get_func(int point) {if (state==HC_COEFF) {fft_coeff_to_func();} return func[point]; }
	double get_coeff(int point) {if (state==HC_FUNC) {fft_func_to_coeff();} return coeff[point]; }

	//operations on the function
	int argmax();
	double valuemax();
	void func_set(int point, double f) {func[point]=f; }
	void func_increment(int point, double f) {func[point]+=f; }
	int normalize(double targetnorm=1.0);
	int reset();
	int scale(double scale);
	int shift(double shift);
	int test();

protected:
	//random number generator
	gsl_rng *rng;
	unsigned int seed;

private:
	// memory management
	bool mem;
	int allocate_mem();
	int free_mem();
};


#define HG_VERBOSE 0
#define HG_LONGTIMEGEN 1000000
#define HG_CONTINUOUS 10000
#define HG_NOTHING 1e-15
#define HG_EXTINCT -9287465
#define HG_BADARG -879564
#define HG_MEMERR -32656845

/**
 * @brief Low-dimensional population evolving on the hypercube_lowd.
 *
 * This class enables simulation of short genomes (\f$L \lesssim 20\f$) but potentially large populations.
 * Random mutation, recombination and selection are supported.
 * A number of properties of the population can be obtained using methods of this class, including:
 * - genotype and allele frequencies;
 * - statistics on fitness and phenotypic traits;
 * - linkage disequilibrium.
 */
class haploid_lowd {
public:
	// public hypercube_lowds
	hypercube_lowd fitness;
	hypercube_lowd population;

	// construction / destruction
	haploid_lowd(int L=1, int rng_seed=0);
	virtual ~haploid_lowd();

	// population parameters (read/write)
	double carrying_capacity;
	double outcrossing_rate;
	bool free_recombination;
	bool circular;				//topology of the chromosome

	// population parameters (read only)
	int L(){return number_of_loci;}
	int get_number_of_loci(){return number_of_loci;}
	double N(){return population_size;}
	double get_population_size(){return population_size;}
	double get_generation(){return long_time_generation+generation;}
	double get_mutation_rate(int locus, int direction) {return mutation_rates[direction][locus];}

	//initialization
	int set_allele_frequencies(double *freq, unsigned long N);
	int set_genotypes(vector <index_value_pair_t> gt);
	int set_wildtype(unsigned long N);

	// modify population
	int set_recombination_rates(double *rec_rates);
	int set_mutation_rate(double m);
	int set_mutation_rate(double m1, double m2);
	int set_mutation_rate(double* m);
	int set_mutation_rate(double** m);

	//evolution
	int evolve(int gen=1);
	int evolve_norec(int gen=1);
	int evolve_deterministic(int gen=1);

	// readout
	// Note: these functions are for the general public and are not expected to be
	// extremely fast. If speed is a major concern, consider subclassing and working
	// with protected methods.

	// genotype readout
	double get_genotype_frequency(int gt){return population.get_func(gt);}
	
	// allele frequencies
	double get_allele_frequency(int locus){return 0.5*(1+(1<<number_of_loci)*population.get_coeff(1<<locus));}
	double get_chi(int locus){return (1<<number_of_loci)*population.get_coeff(1<<locus);}
	double get_moment(int locus1, int locus2){return (1<<number_of_loci)*population.get_coeff((1<<locus1)+(1<<locus2));}
	double get_LD(int locus1, int locus2){return get_moment(locus1, locus2)-get_chi(locus1)*get_chi(locus2);}
	double genotype_entropy();
	double allele_entropy();

	// fitness/phenotype readout
	double get_fitness(int n) {return fitness.get_func(n);}
	stat_t get_fitness_statistics();

	//evolution
	int select();
	int mutate();
	int recombine();
	int resample();

protected:
	//random number generator used for resampling and seeding the hypercube_lowds
	gsl_rng* rng;	//uses the same RNG as defined in hypercube_lowd.h from the  GSL library.
	int seed;	//seed of the rng
	int get_random_seed(); //get random seed from the OS

	//hypercube_lowds that store the distribution of recombinations and the change in the
	//population distribution due to mutations
	hypercube_lowd recombinants;
	hypercube_lowd mutants;
	double** recombination_patters;				// array that holds the probabilities of all possible recombination outcomes for every subset of loci

	// population parameters
	int number_of_loci;
	double population_size;
	int generation;
	double long_time_generation;
	double** mutation_rates;				// the mutation rate can be made locus specific and genotype dependent.

	int calculate_recombinants_free();
	int calculate_recombinants_general();

private:
	// Memory management is private, subclasses must take care only of their own memory
	bool mem;
	int allocate_mem();
	int free_mem();

	// counting reference
	static size_t number_of_instances;
};


class haploid_lowd_test : public haploid_lowd {
public:
	// construction / destruction
	haploid_lowd_test(int L=1, int rngseed=0) : haploid_lowd(L, rngseed){};
	~haploid_lowd_test();

	//testing
	int test_recombinant_distribution();
	int test_recombination(double *rec_rates);
	int mutation_drift_equilibrium(double** mutrates);
};

#endif /* FFPOPSIM_LOWD_H_ */
