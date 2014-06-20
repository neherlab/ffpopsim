/**
 * @file popgen_lowd.h
 * @brief Header file for low-dimensional simulations
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-04-19
 *
 * Copyright (c) 2012-2013, Richard Neher,Fabio Zanini
 * All rights reserved.
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
	bool circular;				//topology of the chromosome

	// population parameters (read only)
	int L(){return number_of_loci;}
	int get_number_of_loci(){return number_of_loci;}
	double N(){return population_size;}
	double get_population_size(){return population_size;}
	double get_generation(){return long_time_generation+generation;}
        void set_generation(double g){if(g > HG_LONGTIMEGEN) {generation = fmod(g, HG_LONGTIMEGEN); long_time_generation = g - generation;} else generation = g;}
	double get_mutation_rate(int locus, int direction) {return mutation_rates[direction][locus];}
	int get_recombination_model(){return recombination_model;}
	double get_recombination_rate(int locus);

	//initialization
	int set_allele_frequencies(double* frequencies, unsigned long N);
	int set_genotypes(vector <index_value_pair_t> gt);
	int set_wildtype(unsigned long N);

	// modify population
	int set_recombination_model(int rec_model);
	int set_recombination_rates(double *rec_rates, int rec_model=-1);
	int set_mutation_rates(double m);
	int set_mutation_rates(double m1, double m2);
	int set_mutation_rates(double* m);
	int set_mutation_rates(double** m);

	//evolution
	int evolve(int gen=1);
	int evolve_norec(int gen=1);
	int evolve_deterministic(int gen=1);

	// readout
	// Note: these functions are for the general public and are not expected to be
	// extremely fast. If speed is a major concern, consider subclassing and working
	// with protected methods.

	// genotype readout
	double get_genotype_frequency(int genotype){return population.get_func(genotype);}
	
	// allele frequencies
	double get_allele_frequency(int locus){return 0.5 * (1 + get_chi(locus));}
	double get_pair_frequency(int locus1, int locus2){return 0.25 * (get_moment(locus1, locus2) - 1) + 0.5 * (get_allele_frequency(locus1) + get_allele_frequency(locus2));}

	double get_chi(int locus){return (1<<number_of_loci)*population.get_coeff(1<<locus);}
	double get_chi2(int locus1, int locus2){return get_moment(locus1, locus2)-get_chi(locus1)*get_chi(locus2);}
	double get_LD(int locus1, int locus2){return 0.25 * get_chi2(locus1, locus2);}
	double get_moment(int locus1, int locus2){return (1<<number_of_loci)*population.get_coeff((1<<locus1)+(1<<locus2));}

        // entropy
	double genotype_entropy();
	double allele_entropy();

	// fitness/phenotype readout
	double get_fitness(int genotype) {return fitness.get_func(genotype);}
	double get_fitness_coefficient(int bitset_loci) {return fitness.get_coeff(bitset_loci);}
	stat_t get_fitness_statistics();

protected:
	//random number generator used for resampling and seeding the hypercube_lowds
	gsl_rng* rng;	//uses the same RNG as defined in hypercube_lowd.h from the  GSL library.
	int seed;	//seed of the rng
	int get_random_seed(); //get random seed from the OS

	//hypercube_lowds that store the distribution of recombinations and the change in the
	//population distribution due to mutations
	hypercube_lowd recombinants;
	hypercube_lowd mutants;
	double** recombination_patterns;	// array that holds the probabilities of all possible recombination outcomes for every subset of loci
	int recombination_model;			//model of recombination to be used

	// population parameters
	int number_of_loci;
	double population_size;
	int generation;
	double long_time_generation;
	double** mutation_rates;		// the mutation rate can be made locus specific and genotype dependent.

	//evolution
	int select();
	int mutate();
	int recombine();
	int resample();

	// recombination
	int set_recombination_rates_general(double *rec_rates);
	int set_recombination_patterns(vector<index_value_pair_t> iv);
	int marginalize_recombination_patterns();
	int set_recombination_rates_single_crossover(double *rec_rates);
	int calculate_recombinants_free();
	int calculate_recombinants_single();
	int calculate_recombinants_general();

private:
	// Memory management is private, subclasses must take care only of their own memory
	bool mem;
	int allocate_mem();
	int free_mem();
	int allocate_recombination_mem(int rec_model);
	int free_recombination_mem();

	// counting reference
	static size_t number_of_instances;
};

#endif /* FFPOPSIM_LOWD_H_ */
