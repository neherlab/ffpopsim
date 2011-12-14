/*
 * popgen.h
 *
 *  Created on: Oct 27, 2010
 *      Author: richard
 */

#ifndef POPGEN_H_
#define POPGEN_H_


#include <iostream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <vector>
#include <bitset>
#include <string>
#include <sstream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <boost/dynamic_bitset.hpp>
using namespace std;

#define MIN(a,b) (a<b)?a:b

#ifndef SAMPLE_H_
#define SAMPLE_H_
#define SAMPLE_ERROR -12312154


class sample {
public:
	int number_of_values;
	double *values;
	double mean;
	double variance;

	int bins;
	bool mem_dis;
	bool mem_values;
	gsl_histogram *distribution;

	bool with_range;
	double range_min;
	double range_max;

	sample();
	virtual ~sample();
	int set_up(int n);
	int set_distribution(int bins=100);
	void set_range(double min, double max) {range_min=min; range_max=max; with_range=true;}
	int calc_mean();
	int calc_variance();
	int calc_distribution();
	int print_distribution(ostream &out);
};

#endif /* SAMPLE_H_ */



#define HC_MEMERR -131545
#define HC_BADARG -131546
#define RNG gsl_rng_taus	//choose the random number generator algorithm, see http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
#define HC_VERBOSE 0			//debugging: if set to one, each function prints out a message into the error stream

#define HC_FUNC 1				//hypercube.func is up-to-date
#define HC_COEFF -1			//hypercube.coeff is up-to-date
#define HC_FUNC_EQ_COEFF 0		//hypercube.func equal hypercube.coeff

using namespace std;

struct index_value_pair{
	int index;
	double val;
};

class hypercube
{
	//dimension of the hypercube
	int dim;

	//random number generator
	gsl_rng *rng;
	unsigned int seed;

	//checks and memory
	bool mem;
	int state;		//takes values HC_FUNC, HC_COEFF, HC_HC_FUNC_EQ_COEFF, depending on the current state of hypercube
	int allocate_mem();
	int free_mem();
	const static double hypercubeversion=0.91;
public:
	double *coeff;		//array holding 2^N coefficients.
						///a entry 0101001101 corresponds to a term with spins at each 1
	double *func;		//array holding the values of the function on the hypercube

	int *order;			//Auxiliary array holding the number of spins, i.e. the number of ones of coeff[k]

	hypercube();
	hypercube(int dim_in, int s=0);
	~hypercube();
	int set_up(int dim_in, int s=0);
	double version(){return hypercubeversion;}

	int gaussian_coefficients(double* vark, bool add=false);
	int additive(double* additive_effects, bool add=false);
	int init_rand_gauss(double sigma, bool add=false);
	int init_list(vector<index_value_pair> iv, bool add=false);
	int init_coeff_list(vector <index_value_pair> iv, bool add=false);
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
	int get_state() {return state;};
	unsigned int get_seed() {return seed;};
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
};


#define HG_LONGTIMEGEN 1000000
#define HG_CONTINUOUS 10000
#define HG_NOTHING 1e-15
#define HG_EXTINCT -9287465
#define HG_MEMERR -32656845

class haploid_gt_dis
{
	//hypercubes that store the distribution of recombinations and the change in the
	//population distribution due to muations
	hypercube recombinants;
	hypercube mutants;
	//array that holds the probabilities of all possible recombination outcomes for every subset of loci
	double** recombination_patters;

	double population_size;
	int number_of_loci;
	int generation;
	double long_time_generation;
	//mutation rate could be made locus specific and genotype dependent.
	double** mutation_rates;
	double outcrossing_rate;
	bool free_recombination;

	//random number generator used for resampling and seeding the hypercubes
	gsl_rng* rng;	//uses the same RNG as defined in hypercube.h from the  GSL library.
	int seed;	//seed of the rng

	int allocate_mem();
	int free_mem();
	bool mem;
	const static double haploidversion=0.61;
public:
	hypercube fitness;
	hypercube population;
	//setting up
	haploid_gt_dis();
	~haploid_gt_dis();
	haploid_gt_dis(int nloci, double popsize, int rngseed=0);
	int setup(int nloci, double popsize, int rngseed=0);
	double version(){return haploidversion;}

	//initialization
	int init_frequencies(double *freq);
	int init_genotypes(vector <index_value_pair> gt);
	int set_recombination_rates(double *rec_rates);
	void set_population_size(double popsize){population_size=popsize;}
	void set_outcrossing_rate(double orate){outcrossing_rate=orate;}
	int set_mutation_rate(double m);
	int set_mutation_rate(double m1, double m2);
	int set_mutation_rate(double* m);
	int set_mutation_rate(double** m);

	//evolution
	int evolve(int gen=1);
	int evolve_norec(int gen=1);
	int evolve_deterministic(int gen=1);
	int select();
	int mutate();
	int recombine();
	int resample(double n=0.0);
	int calculate_recombinants();
	int calculate_recombinants_free();
	int calculate_recombinants_general();

	//analyze and access population
	double get_genotype_frequency(int gt){return population.get_func(gt);}
	double get_allele_frequency(int locus){return 0.5*(1+(1<<number_of_loci)*population.get_coeff(1<<locus));}
	double get_chi(int locus){return (1<<number_of_loci)*population.get_coeff(1<<locus);}
	double get_moment(int locus1, int locus2){return (1<<number_of_loci)*population.get_coeff((1<<locus1)+(1<<locus2));}
	double get_LD(int locus1, int locus2){return get_moment(locus1, locus2)-get_chi(locus1)*get_chi(locus2);}
	double get_generation(){return long_time_generation+generation;}
	double genotype_entropy();
	double allele_entropy();
	double fitness_mean();
	double fitness_variance();
	double get_population_size(){return population_size;}
	double N(){return population_size;}
	int L(){return number_of_loci;}
	double get_mutation_rate(int locus, int direction) {return mutation_rates[direction][locus];}
	double get_outcrossing_rate() {return outcrossing_rate;}
	//testing
	int test_recombinant_distribution();
	int test_recombination(double *rec_rates);
	int mutation_drift_equilibrium(double** mutrates);
};


#define HC_MEMERR -131545
#define HC_BADARG -131546
#define HC_VERBOSE 0


using namespace std;

class coeff
{
public:
	int order;
	double value;
	int *loci;
	int *words;
	int *bits;
	coeff(double value_in, vector <int> loci_in){
		value=value_in;
		order=loci_in.size();
		loci=new int [order];
		words=new int [order];
		bits=new int [order];
		for (int i=0; i<order; i++) loci[i]=loci_in[i];
	}
	~coeff(){
		//delete [] loci;
		//delete [] words;
		//delete [] bits;
	}
};

class coeff_single_locus{
public:
	double value;
	int locus;
	int word;
	int bit;
	coeff_single_locus(double value_in, int locus_in){
		value=value_in;
		locus=locus_in;
	}

};

class hypercube_function
{
public:
	int dim;
	double hypercube_mean;
	vector <coeff_single_locus> coefficients_single_locus;
	vector <coeff> coefficients_epistasis;
	double epistatic_std;
	int rng_offset;

	//random number generator
	gsl_rng *rng;
	unsigned int seed;

	bool hcube_allocated;
	bool mem;
	int allocate_mem();
	int free_mem();

	hypercube_function();
	~hypercube_function();
	int set_up(int dim_in,  int s=0);

	unsigned int get_seed() {return seed;};
	double get_func(boost::dynamic_bitset<> *gt);
	double get_func_words(int *gt, int n_o_w);
	double get_additive_coefficient(int locus);
	int add_coefficient(double value, vector <int> loci);
	int add_coefficient(double value, vector <int> loci, int n_o_w, int* n_o_b);
	int set_additive_coefficient(double value, int locus, int expected_locus=-1);
	int set_random_epistasis_strength(double sigma);
};


#ifndef HAPLOIDPOPULATION_H_
#define HAPLOIDPOPULATION_H_
#define HP_VERBOSE 0
#define NO_GENOTYPE -1
#define HP_BADARG -879564
#define HP_MEMERR -986465
#define HP_MINAF 0.02
#define HP_NOTHING 1e-12
#define HP_RANDOM_SAMPLE_FRAC 0.01
#define FREE_RECOMBINATION 1
#define CROSSOVERS 2
#define WORD_LENGTH 20

/* Use dynamic bitsets because they require little memory */
struct gt {
	boost::dynamic_bitset<> genotype;
	vector <double> trait;
	double fitness;
	int clone_size;
	gt(int n_traits){trait.resize(n_traits);}
};

struct stat{
	double mean;
	double variance;
};


struct label_gt_pair{ int gt; double label; };


class haploid_clone {
	int number_of_traits;
	int number_of_individuals;		//maximal population size in terms of memory allocated to hold genotypes
	int pop_size;				//actual population size
	int target_pop_size;			//average population size
	int scratch;				//variable by how much the memory for offsprings exceeds the number_of_individuals (1+scratch)*..
	int generation;

	double mutation_rate;			//rate of mutation per locus per generation
	double outcrossing_probability;
	double crossover_rate;
	bool circular;				//topology of the chromosome
	int recombination_model;		//model of recombination to be used

	int number_of_loci;			//total number of loci
	vector <int> sex_gametes;		//array holding the indices of gametes
	vector <int> random_sample;		//array holding the indices of gametes

	int *genome;				//Auxiliary array holding the positions along the genome
	int *crossovers;

	stat fitness_stat;			//structure holding the fitness values of the individuals in the population
	stat *trait_stat;
	double **trait_covariance;

	double *allele_frequencies;
	double *gamete_allele_frequencies;
	double *chi1;				//symmetric allele frequencies
	double **chi2;				//symmetric two locus correlations

	bool mem;
	bool cumulants_mem;
	bool evolve_rec_rates;

	int allocate_mem();
	int free_mem();
	gsl_rng* evo_generator;
	gsl_rng* label_generator;
	int seed;

	int recombine_crossover(int parent1, int parent2, int ng);
	boost::dynamic_bitset<> reassortment_pattern();
	boost::dynamic_bitset<> crossover_pattern();
	int recombine(int parent1, int parent2);
	double chemical_potential();

	void flip_single_locus(int individual, int locus);

public:
	hypercube_function *trait;		// genotype to fitness map
	vector <gt> *current_pop;		// genotypes is a [number_of_individuals*number_of_words] array storing the alleles of each individual
	vector <gt> *new_pop;			// newgenotypes is a [number_of_individuals*number_of_words] array storing the alleles of each individual
	vector <gt> pop2;			// genotypes is a [number_of_individuals*number_of_words] array storing the alleles of each individual
	vector <gt> pop1;			// newgenotypes is a [number_of_individuals*number_of_words] array storing the alleles of each individual

	haploid_clone();
	virtual ~haploid_clone();
	int set_up(int N_in,int L,  int rng_seed=0, int number_of_traits=1);
	void set_target_pop_size(int tgps){target_pop_size=tgps;}
	void set_mutation_rate(double mu){mutation_rate=mu;}
	void set_outcrossing_probability(double r){outcrossing_probability=r;}
	void set_fixed_rec_rate(double r){outcrossing_probability=r;}		// deprecated
	void set_crossover_rate(double c){crossover_rate=c;}
	void set_recombination_model(int c) {recombination_model=c;}
	void set_circular(bool c) {circular=c;}
	void produce_random_sample(int size);
	int init_genotypes(double *nu, int n_o_genotypes=0);
	int init_genotypes_diverse(int n_o_genotypes, int no_copies=1);
	int init_genotypes(int n_o_genotypes=-1);

	int evolve();
	int bottleneck(int size_of_bottleneck);
	void mutate();
	int select_gametes();
	int add_recombinants();
	int swap_populations();
	int new_generation();
	int flip_single_locus(int locus);
	void shuffle_genotypes();
	void update_fitness();
	void calc_stat();
	void calc_fit();
	void calc_fitness_stat();
	void calc_trait_stat();
	void calc_individual_fitness(gt *tempgt);
	void calc_everybodies_traits();
	virtual void calc_fitness_from_traits(gt *tempgt){tempgt->fitness = tempgt->trait[0];}
	double fitness_mean() {return fitness_stat.mean;}
	double fitness_var() {return fitness_stat.variance;}
	double trait_mean(int t=0) {return trait_stat[t].mean;}
	double trait_var(int t=0) {return trait_stat[t].variance;}
	double trait_cov(int t1, int t2) {return trait_covariance[t1][t2];}
	double get_multi_point_frequency(vector <int> loci);
	double get_pair_frequency(int locus1, int locus2);
	vector <double> get_pair_frequencies(vector < vector <int> > *loci);
	double get_allele_frequency(int l) {return allele_frequencies[l];}
	double get_chi(int l) {return 2*allele_frequencies[l]-1.0;}
	double get_number_of_clones(){return current_pop->size();}
	int N() {return pop_size;}
	double mu() {return mutation_rate;}
	int get_pop_size() {return pop_size;}		//Deprecated
	int get_target_pop_size() {return target_pop_size;}
	int L(){return number_of_loci;}

	double diversity_mean();
	double divergence_mean();
	unsigned int distance_Hamming(boost::dynamic_bitset<> genotype, boost::dynamic_bitset<> genotype1);
	vector <unsigned long> partition_cumulative();

	//int get_genotype(int i) {return genotypes[i].genotype;}
	string get_genotype_string(int i);
	void add_genotypes(boost::dynamic_bitset<> newgt,  int n);
	int add_fitness_coefficient(double value, vector <int> loci, int traitnumber=0){return trait[traitnumber].add_coefficient(value, loci);}
	void clear_fitness_function(){for (int t=0; t<number_of_traits; t++){trait[t].coefficients_single_locus.clear(); trait[t].coefficients_epistasis.clear();}}
	int random_clone();

	double get_fitness(int n) {return (*current_pop)[n].fitness;}
	double get_trait(int n, int t) {return (*current_pop)[n].trait[t];}
	double get_max_fitness();
	int get_generation(){return generation;}
	int print_allele_frequencies(ostream &out);
	int read_ms_sample(istream &gts, int skip_locus, int multiplicity);
	int read_ms_sample_sparse(istream &gts, int skip_locus, int multiplicity, int distance);
};

class haploid_population {
	int number_of_individuals;	//maximal population size in terms of memory allocated to hold genotypes
	int pop_size;				//actual population size
	int target_pop_size;		//average population size
	int scratch;				//variable by how much the memory for offsprings exceeds the number_of_individuals (1+scratch)*..
	int generation;

	double mutation_rate; 		//rate of mutation per locus per generation
	double fixed_outcrossing_rate;		//fixed recombination rate, in case no recombination modifiers needed
	bool circular;				//topology of the chromosome
	int recombination_model;	//model of recombination to be used

	int number_of_loci;			//total number of loci
	int number_of_words;		//number of integers used to store one genotype
	int *number_of_bits;		//number of bits used in each word
	int *genotypes;				//genotypes is a [number_of_individuals*number_of_words] array storing the alleles of each individual
	int *new_genotypes;			//newgenotypes is a [number_of_individuals*number_of_words] array storing the alleles of each individual
	int *sex_gametes;				//array holding the indices of gametes
	int *asex_gametes;				//array holding the indices of gametes
	int n_asex_gametes;
	int n_sex_gametes;

	int *clone_size_distribution;
	int number_of_clones;

	int *genome;				//Auxiliary array holding the positions along the genome
	int *crossovers;			//

	sample *fitness_values;		//structure holding the fitness values of the individuals in the population
	sample *new_fitness_values;	//same as above, used for swapping
	sample *outcrossing_rates;			//structure holding the recombination rate of the individuals in the population
	sample *new_outcrossing_rates;		//same as above, used for swapping
	sample *gt_labels;			//structure holding the labels of the genotype, this is used for genotype comparison and genotype sorting
	sample *new_gt_labels;		//same as above, used for swapping

	double *allele_frequencies;
	double *gamete_allele_frequencies;
	double *chi1;				//symmetric allele frequencies
	double **chi2;				//symmetric two locus correlations

	bool mem;
	bool clone_size_mem;
	bool cumulants_mem;
	bool evolve_outcrossing_rates;

	int allocate_mem();
	int free_mem();
	gsl_rng* evo_generator;
	gsl_rng* label_generator;
	int seed;
	int* numbers;

	int recombine_crossover(int parent1, int parent2, int ng);
	int recombine_free(int parent1, int parent2, int ng);
	double chemical_potential();

	int index(int individual, int word) {return individual*number_of_words+word;}
	int locus_word(int locus) {return locus / WORD_LENGTH;}
	int locus_bit(int locus) {return locus % WORD_LENGTH;}
	void flip_single_locus(int individual, int locus);

public:
	hypercube_function fitness;	//genotype to fitness map
	hypercube_function outcrossing_rate;	//genotype to recombination map


	haploid_population();
	virtual ~haploid_population();
	int set_up(int N_in,int L,  int rng_seed=0);
	void set_target_pop_size(int tgps){target_pop_size=tgps;}
	void set_mutation_rate(double mu){mutation_rate=mu;}
	void set_fixed_outcrossing_rate(double r){fixed_outcrossing_rate=r; evolve_outcrossing_rates=false;}
	void set_fixed_rec_rate(double r){fixed_outcrossing_rate=r; evolve_outcrossing_rates=false;} //compatibility function kept after making rec_rate outcrossing_rate
	void set_recombination_model(int c) {recombination_model=c;}
	void set_circular(bool c) {circular=c;}
	int allocate_clone_size_distribution();
	int set_evolve_outcrossing_rates();

	int init_genotypes(double *nu, int n_o_genotypes=0);
	int init_genotypes_diverse(int n_o_genotypes, int no_copies=1);
	int init_genotypes(int n_o_genotypes);

	int evolve();
	int bottleneck(int size_of_bottleneck);
	void mutate();
	int select_gametes();
	int new_generation();
	int flip_single_locus(int locus);
	void shuffle_genotypes();

	void calc_stat();
	void calc_rec();
	void calc_fit();
	void calc_gt_labels();
	void calc_individual_fitness(int individual);
	double calc_label(int* individual);
	int calc_clone_size_distribution();

	double fitness_mean() {return fitness_values->mean;}
	double fitness_var() {return fitness_values->variance;}
	double rec_mean() {return outcrossing_rates->mean;}
	double rec_var() {return outcrossing_rates->variance;}
	int get_pop_size() {return pop_size;}
	double get_multi_point_frequency(vector <int> loci);
	double get_allele_frequency(int l) {return allele_frequencies[l];}
	double get_chi(int l) {return 2*allele_frequencies[l]-1.0;}
	double get_number_of_clones(){if (clone_size_mem) return number_of_clones; else return HP_BADARG;}
	double get_clone_size(int c){if (clone_size_mem and c<number_of_clones) return clone_size_distribution[2*c]; else return HP_BADARG;}
	double get_clone_individual(int c){if (clone_size_mem and c<number_of_clones) return clone_size_distribution[2*c+1]; else return HP_BADARG;}
	int L(){return number_of_loci;}
	int N(){return pop_size;}
	double mu(){return mutation_rate;}

	int get_genotype(int i, int w) {if (i<pop_size and w<number_of_words) return genotypes[index(i,w)]; else return NO_GENOTYPE;}
	string get_genotype_string(int i);
	int remove_last_genotypes(int m);
	int add_genotype(int *gt);
	int add_genotypes(int *gt, int n);
	int add_fitness_coefficient(double value, vector <int> loci){return fitness.add_coefficient(value, loci, number_of_words, number_of_bits);}
	int add_rec_coefficient(double value, vector <int> loci){return outcrossing_rate.add_coefficient(value, loci, number_of_words, number_of_bits);}
	void clear_fitness_function(){fitness.coefficients_single_locus.clear(); fitness.coefficients_epistasis.clear();}

	double get_fitness(int n) {return fitness_values->values[n];}
	double get_rec(int n) {return outcrossing_rates->values[n];}
	double get_label(int n) {return gt_labels->values[n];}
	double get_max_fitness();
	int get_generation(){return generation;}
	int print_allele_frequencies(ostream &out);
	int read_ms_sample(istream &gts, int initial_locus, int multiplicity=1);
};


#endif /* HAPLOIDPOPULATION_H_ */


#endif /* POPGEN_H_ */

