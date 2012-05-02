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

/*
%apply (int DIM1, double* INPLACE_ARRAY1) {(int gts2, double* frequencies)};
%apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int gts, int loci, double* hoprates)};
%apply (int DIM1, double* IN_ARRAY1) {(int loci, double* fitness_additive)};
*/

/* Specify external objects */
struct stat_t {
	double mean;
	double variance;
};

/*TODO: add C++ wrappers for the pointers (set/get methods) */
class hypercube_function {
public:
	int dim;
	double hypercube_mean;
	vector <coeff_single_locus_t> coefficients_single_locus;
	vector <coeff_t> coefficients_epistasis;
	double epistatic_std;
	int rng_offset;

	// setting up
	hypercube_function();
	hypercube_function(int dim_in, int s=0);
	virtual ~hypercube_function();
	int set_up(int dim_in,  int s=0);

	// methods
	unsigned int get_seed() {return seed;};
	double get_func(boost::dynamic_bitset<> *genotype);
	double get_additive_coefficient(int locus);
	int set_additive_coefficient(double value, int locus, int expected_locus=-1);
	int add_coefficient(double value, vector <int> loci);
	int set_random_epistasis_strength(double sigma);
};


/*TODO: add C++ wrappers for the pointers (set/get methods) */
struct clone_t {
	boost::dynamic_bitset<> genotype;
	vector <double> trait;
	double fitness;
	int clone_size;
	clone_t(int n_traits){trait.resize(n_traits);}
};

class haploid_clone {
public:
	// population parameters (read only)
	int get_generation(){return generation;}
	int get_number_of_loci(){return number_of_loci;}
	int get_pop_size() {return pop_size;}
	int get_number_of_clones(){return current_pop->size();}

        //evolve
	int bottleneck(int size_of_bottleneck);

	// population parameters (read/write)
	int target_pop_size;			// target (average) population size
	double mutation_rate;			// rate of mutation per locus per generation
	double outcrossing_probability;		// probability of having sex
	double crossover_rate;			// rate of crossover during sex
	int recombination_model;		//model of recombination to be used
	bool circular;				//topology of the chromosome

	// allele frequencies
	double get_allele_frequency(int l) {return allele_frequencies[l];}
	double get_pair_frequency(int locus1, int locus2);

	// fitness/phenotype readout
	double get_fitness(int n) {calc_individual_fitness(&((*current_pop)[n])); return (*current_pop)[n].fitness;}
	double get_trait(int n, int t=0) {calc_individual_traits(&((*current_pop)[n])); return (*current_pop)[n].trait[t];}

};



%include "hivpython.h"

