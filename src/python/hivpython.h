/**
 * @file hivpopulation.h
 * @brief Header file for a typical HIV population (subclass of haploid_clone)
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-23
 */
#ifndef HIVPYTHON_H_
#define HIVPYTHON_H_

#include "hivpopulation.h"

#define HIVPYTHON_VERBOSE 0

/**
 * @brief HIV population with facultative drug treatment (Python2 bindings)
 *
 * This class is an interface class for the Python bindings (like a dressing screen).
 * It contains an hivpopulation instance as a private attribute, so that problematic
 * bindings are avoided, and exposes only a few, Python-friendly methods.
 *
 * Note: whenever a baseclass method is used "as is", it is declared as public via C++'s
 * "using" mechanism. The corresponsing line of the baseclass header file is included by
 * hand in the SWIG interface file.
 */
class hivpython: private hivpopulation {
private:

public:
	// constructors/destructors
	hivpython();
	~hivpython();

	// we mask some baseclass methods by equivalent ones that are described in the .cpp file
	// this way we can selectively include code into SWIG
	int set_up(int N_in, int rng_seed=0, double mutrate=3e-5, double coinfection_rate=1e-2, double crossover_rate=1e-3){return hivpopulation::set_up(N_in, rng_seed, mutrate, coinfection_rate, crossover_rate);}
	int init_genotypes(int n_o_genotypes=-1){return haploid_clone::init_genotypes(n_o_genotypes);}
	int init_genotypes(double IN_ARRAY1[HIVGENOME], int n_o_genotypes=0);


	// population parameters (read only)
	using haploid_clone::get_generation;
	using haploid_clone::get_number_of_loci;
	using haploid_clone::get_pop_size;
	using haploid_clone::get_number_of_clones;

	// population parameters (read/write)
	using haploid_clone::target_pop_size;
	using haploid_clone::mutation_rate;
	using haploid_clone::outcrossing_probability;
	using haploid_clone::crossover_rate;
	using haploid_clone::recombination_model;
	using haploid_clone::circular;

	// evolve
	int evolve(int gen=1);
	using haploid_clone::bottleneck;

	// genotype readout
	// Note: one cannot draw several clones in one shot because of limitations in
	// the SWIG Numpy interface (might be possible in the future). This is still exendible
	// on the python side if wished (but might be expensive).
	using haploid_clone::random_clone;
	void random_clones(int DIM1, unsigned int * ARGOUT_ARRAY1);
	void get_genotype(unsigned int i, unsigned short ARGOUT_ARRAY1[HIVGENOME]);
	int distance_Hamming(unsigned int clone1, unsigned int clone2){return haploid_clone::distance_Hamming(clone1, clone2);}
	int distance_Hamming(unsigned int clone1, unsigned int clone2, int DIM1, int DIM2, unsigned int * IN_ARRAY2, int every=1);
	stat_t get_diversity_statistics(unsigned int n_sample=1000){haploid_clone::calc_stat(); return haploid_clone::get_diversity_statistics(n_sample);}
	stat_t get_divergence_statistics(unsigned int n_sample=1000){haploid_clone::calc_stat(); return haploid_clone::get_divergence_statistics(n_sample);}
//FIXME	int get_divergence_histogram(int DIM1, unsigned int * INPLACE_ARRAY1, int DIM1, unsigned int * INPLACE_ARRAY1, int DIM1, int DIM2, unsigned int * IN_ARRAY2, unsigned int every=1, unsigned int n_sample=1000);

	// fitness/phenotype readout
	using haploid_clone::get_fitness;
	using haploid_clone::get_trait;
	void get_fitnesses(int DIM1, double* ARGOUT_ARRAY1);

	// allele frequencies
	using haploid_clone::get_allele_frequency;
	using haploid_clone::get_pair_frequency;
	void get_allele_frequencies(double ARGOUT_ARRAY1[HIVGENOME]);

	// treatment (set/get)
	void set_treatment(double t){hivpopulation::set_treatment(t);}
	double get_treatment(){return hivpopulation::get_treatment();}
	void calc_fitness_from_traits(int l){hivpopulation::calc_fitness_from_traits(&((*current_pop)[l]));}

	// stream I/O
	int read_selection_coefficients(char *model);
	int read_resistance_coefficients(char *model);
	int write_genotypes(char * filename, int sample_size, char * gt_label=NULL, int start=0, int length=0);
};

#endif /* HIVPYTHON_H_ */
