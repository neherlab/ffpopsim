/**
 * @file popgen_highd.h
 * @brief Header file for high-dimensional simulations
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-04-19
 *
 * Copyright (c) 2012, Richard Neher, Fabio Zanini
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
 *
 * HP_VERBOSE: degree of verbosity of haploid_highd. Levels:
 * - 0: no messages
 * - 1: most messages (enter/exit function)
 * - 2: all messages
 */
#ifndef FFPOPSIM_HIGHD_H_
#define FFPOPSIM_HIGHD_H_
#include "ffpopsim_generic.h"

#define HCF_MEMERR -131545
#define HCF_BADARG -131546
#define HCF_VERBOSE 0
#define WORDLENGTH 28 	//length used to chop bitsets into words

using namespace std;


/**
 * @brief Trait coefficient for a set of loci.
 *
 * This struct is used in the hypercube_highd class for saving trait coefficients. See also hypercube_highd::add_coefficient.
 *
 * Note: words and bits are deprecated.
 */
struct coeff_t {
	int order;
	double value;
	int *loci;
	coeff_t(double value_in, vector <int> loci_in){
		value=value_in;
		order=loci_in.size();
		loci=new int [order];
		for (int i=0; i<order; i++) loci[i]=loci_in[i];
	}
};

/**
 * @brief Trait coefficient for a single locus.
 *
 * Note: word and bit are deprecated. For this reason, this class is actually equivalent
 * to the index_value_pair_t struct, with a mandatory constructor.
 */
struct coeff_single_locus_t {
	double value;
	int locus;
	coeff_single_locus_t(double value_in, int locus_in) : value(value_in), locus(locus_in) {};
};

/**
 * @brief Hypercube class for high-dimensional simulations.
 *
 * This class is used for representing properties of genotypes.
 * Unlike the low-dimensional sister class, no attempt is made to monitor the whole hypercube,
 * because the dimension the space of genotypes of length L is \f$2^L\f$.
 * As a consequence, no Fourier transformations are implemented.
 *
 * The main use of this class is either for pointwise assigment of values
 * (for instance, we assign fitness only for the actually observed genotypes)
 * or for low-order Fourier coefficients, like like main fitness effects and two-site epistasis,
 * which scale like L and \f$L^2\f$, respectively.
 *
 */
class hypercube_highd {
private:
	//random number generator
	gsl_rng *rng;
	unsigned int seed;

	// memory management
	bool hcube_allocated;
	bool mem;
	int allocate_mem();
	int free_mem();

public:
	int dim;
	double hypercube_mean;
	vector <coeff_single_locus_t> coefficients_single_locus;
	vector <coeff_t> coefficients_epistasis;
	double epistatic_std;
	int rng_offset;

	// setting up
	hypercube_highd();
	hypercube_highd(int dim_in, int s=0);
	virtual ~hypercube_highd();
	int set_up(int dim_in,  int s=0);

	// get methods
	unsigned int get_dim(){return dim;}
	unsigned int get_seed() {return seed;};
	double get_func(boost::dynamic_bitset<> *genotype);
	double get_additive_coefficient(int locus);

	// change the hypercube
	void reset();
	int set_additive_coefficient(double value, int locus, int expected_locus=-1);
	int add_coefficient(double value, vector <int> loci);
	int set_random_epistasis_strength(double sigma);
};


// Control constants
#define HP_VERBOSE 0
#define NO_GENOTYPE -1
#define HP_MINAF 0.02
#define MAX_DELTAFITNESS 8
#define MAX_POPSIZE 500000
#define HP_NOTHING 1e-12
#define HP_RANDOM_SAMPLE_FRAC 0.01
#define FREE_RECOMBINATION 1
#define CROSSOVERS 2

// Error Codes
#define HP_BADARG -879564
#define HP_MEMERR -986465
#define HP_EXPLOSIONWARN 4
#define HP_EXTINCTERR 5
#define HP_NOBINSERR 6
#define HP_WRONGBINSERR 7


/**
 * @brief clone with a single genotype and a vector of phenotypic traits.
 *
 * Note: it uses dynamic bitsets because they require little memory.
 */
struct clone_t {
	boost::dynamic_bitset<> genotype;
	vector <double> trait;
	double fitness;
	int clone_size;
	clone_t(int n_traits=0) : genotype(boost::dynamic_bitset<>(0)), trait(n_traits, 0), fitness(0), clone_size(0) {};

        // Comparison operators check fitness first, genome (big endian) last
	bool operator==(const clone_t &other) const {return (fitness == other.fitness) && (genotype == other.genotype);}
        bool operator!=(const clone_t &other) const {return (fitness != other.fitness) || (genotype != other.genotype);}
	bool operator<(const clone_t &other) const {
                if(fitness < other.fitness) return true;
                else if (fitness > other.fitness) return false;
                else {
                        for(size_t i=0; i < genotype.size(); i++) {
                                if((!genotype[i]) && (other.genotype[i])) return true;
                                else if((genotype[i]) && (!other.genotype[i])) return false;
                        }
                        return false;
                }
        }
	bool operator>(const clone_t &other) const {
                if(fitness > other.fitness) return true;
                else if (fitness < other.fitness) return false;
                else {
                        for(size_t i=0; i < genotype.size(); i++) {
                                if((genotype[i]) && (!other.genotype[i])) return true;
                                else if((!genotype[i]) && (other.genotype[i])) return false;
                        }
                        return false;
                }
        }
};

/**
 * @brief Population class for high-dimensional simulations.
 *
 * This class is the main object storing the state of and enabling the manipulation of populations with long genomes (\f$L \gtrsim 20\f$).
 *
 * Both asexual and sexual populations can be simulated. Since asexual populations under selection are often structured as a small list of
 * large clones, the class stores the clones and their sizes instead of the individuals.
 *
 * Class methods give access to a variety of information on the population, including:
 * - the fitness distribution;
 * - summary statistics of fitness and other phenotypic trits;
 * - genetic structure (linkage disequilibrium, allele frequencies, number of clones).
 */
class haploid_highd {
public:
	// genotype to traits maps, which in turn are used in the trait-to-fitness map
	hypercube_highd *trait;

	// pointers to the current population and the temporary one
	vector <clone_t> *current_pop;
	vector <clone_t> *new_pop;

	// construction / destruction
	haploid_highd(int L=0, int rng_seed=0, int number_of_traits=1);
	virtual ~haploid_highd();

	// population parameters (read/write)
	int carrying_capacity;			// carrying capacity of the environment (pop size)
	double mutation_rate;			// rate of mutation per locus per generation
	double outcrossing_rate;		// probability of having sex
	double crossover_rate;			// rate of crossover during sex
	int recombination_model;		//model of recombination to be used
	bool circular;				//topology of the chromosome

	// population parameters (read only)
	int L(){return number_of_loci;}
	int get_number_of_loci(){return number_of_loci;}
	int N(){return population_size;}
	int get_population_size() {return population_size;}
	int get_generation(){return generation;}
	int get_number_of_clones(){return current_pop->size();}
	int get_number_of_traits(){return number_of_traits;}

	// initialization
	int set_allele_frequencies(double *freq, unsigned long N);
	int set_genotypes(vector <genotype_value_pair_t> gt);
	int set_wildtype(unsigned long N);

	// modify population
	void add_genotypes(boost::dynamic_bitset<> newgt, int n);
	void flip_single_locus(unsigned int clonenum, int locus);

	// modify traits
	int add_trait_coefficient(double value, vector <int> loci, int traitnumber=0){return trait[traitnumber].add_coefficient(value, loci);}
	void clear_traits(){for(int t=0; t<number_of_traits; t++){trait[t].reset();}}
	void set_random_trait_epistasis(double epistasis_std,int traitnumber=0){trait[traitnumber].epistatic_std=epistasis_std;}

	// modify fitness (shortcuts: they only make sense if number_of_traits=1)
	int add_fitness_coefficient(double value, vector <int> loci){if(number_of_traits>1) return HP_BADARG; return add_trait_coefficient(value, loci, 0);}
	void clear_fitness(){if(number_of_traits>1){if(HP_VERBOSE) cerr<<"What do you mean by fitness?"<<endl; throw (int)HP_BADARG;} clear_traits();}
	void set_random_epistasis(double epistasis_std){trait[0].epistatic_std=epistasis_std;}

	// evolution
	int evolve(int gen=1);	
	int bottleneck(int size_of_bottleneck);

	// update traits and fitness and calculate statistics
	void calc_stat();
	void unique_clones();
	

	// readout
	// Note: these functions are for the general public and are not expected to be
	// extremely fast. If speed is a major concern, consider subclassing and working
	// with protected methods.

	// random clones
	int random_clone();
	int random_clones(unsigned int n_o_individuals, vector <int> *sample);

	// genotype readout
	string get_genotype_string(unsigned int i){string gts; boost::to_string((*current_pop)[i].genotype, gts); return gts;}
	int distance_Hamming(unsigned int clone1, unsigned int clone2, vector <unsigned int *> *chunks=NULL, unsigned int every=1){return distance_Hamming((*current_pop)[clone1].genotype, (*current_pop)[clone2].genotype, chunks, every);}
	int distance_Hamming(boost::dynamic_bitset<> gt1, boost::dynamic_bitset<> gt2, vector<unsigned int *> *chunks=NULL, unsigned int every=1);
	stat_t get_diversity_statistics(unsigned int n_sample=1000);
	stat_t get_divergence_statistics(unsigned int n_sample=1000);

	// allele frequencies
	double get_allele_frequency(int l) {return allele_frequencies[l];}
	double get_pair_frequency(int locus1, int locus2);
	vector <double> get_pair_frequencies(vector < vector <int> > *loci);
	double get_chi(int l) {return 2*allele_frequencies[l]-1.0;}
	double get_participation_ratio(){return participation_ratio;};
//	double get_multi_point_frequency(vector <int> loci);

	// fitness/phenotype readout
	double get_fitness(int n) {calc_individual_fitness(&((*current_pop)[n])); return (*current_pop)[n].fitness;}
	double get_trait(int n, int t=0) {calc_individual_traits(&((*current_pop)[n])); return (*current_pop)[n].trait[t];}
	stat_t get_fitness_statistics(){update_fitness(); calc_fitness_stat(); return fitness_stat;}
	stat_t get_trait_statistics(int t=0){calc_trait_stat(); return trait_stat[t];}
	double get_trait_covariance(int t1, int t2) {calc_trait_stat(); return trait_covariance[t1][t2];}
	double get_max_fitness();

	// histograms
	int get_divergence_histogram(gsl_histogram **hist, unsigned int bins=10, vector <unsigned int *> *chunks=NULL, unsigned int every=1, unsigned int n_sample=1000);
	int get_diversity_histogram(gsl_histogram **hist, unsigned int bins=10, vector <unsigned int *> *chunks=NULL, unsigned int every=1, unsigned int n_sample=1000);
	int get_fitness_histogram(gsl_histogram **hist, unsigned int bins=10, unsigned int n_sample=1000);

	// stream I/O
	int print_allele_frequencies(ostream &out);
	int read_ms_sample(istream &gts, int skip_locus, int multiplicity);
	int read_ms_sample_sparse(istream &gts, int skip_locus, int multiplicity, int distance);

protected:
	// random number generator
	gsl_rng* evo_generator;
	gsl_rng* label_generator;
	int seed;
	int get_random_seed();
	vector <int> random_sample;
	void produce_random_sample(int size=1000);

	// population parameters
	int number_of_loci;
	int population_size;
	int number_of_traits;
	int scratch;			//variable by how much the memory for offsprings exceeds the number_of_individuals (1+scratch)*..
	int generation;

	// evolution
	int mutate();
	int select_gametes();
	double relaxation_value();
	double get_logmean_expfitness(double fitness_max);	// Log of the population exp-average of the fitness: log[<exp(F)>_{population}]
	
	int flip_single_locus(int locus);
	void shuffle_genotypes();
	int swap_populations();
	int new_generation();

	// clone structure
	int partition_cumulative(vector <unsigned int> &partition_cum);

	// allele_frequencies
	double *allele_frequencies;
	double *gamete_allele_frequencies;
	double *chi1;				//symmetric allele frequencies
	double **chi2;				//symmetric two locus correlations
	void calc_allele_freqs();

	// recombination details
	double outcrossing_rate_effective;
	int *genome;				//Auxiliary array holding the positions along the genome
	int *crossovers;
	boost::dynamic_bitset<> reassortment_pattern();
	boost::dynamic_bitset<> crossover_pattern();
	vector <int> sex_gametes;		//array holding the indices of gametes
	int add_recombinants();
	int recombine(int parent1, int parent2);
	int recombine_crossover(int parent1, int parent2, int ng);

	// fitness and traits
	double participation_ratio;
	stat_t fitness_stat;
	stat_t *trait_stat;
	double **trait_covariance;
	void update_traits();
	void update_fitness();
	void calc_fitness_stat();
	void calc_trait_stat();
	void calc_individual_traits(clone_t *tempgt);
	void calc_individual_fitness(clone_t *tempgt);
	virtual void calc_individual_fitness_from_traits(clone_t *tempgt){tempgt->fitness = tempgt->trait[0];}	// this must be virtual, because the fitness landscape on the (genotype x phenotype) space can be wild (here fitness IS the only trait)

private:
	// Memory management is private, subclasses must take care only of their own memory
	bool mem;
	bool cumulants_mem;
	int allocate_mem();
	int free_mem();

	// These two vectors are referenced by subclasses and programs using their pointers,
	// current_pop and new_pop, which are public. In fact, these vectors are created only
	// to ensure their memory is released upon destruction of the class.
	vector <clone_t> current_pop_vector;
	vector <clone_t> new_pop_vector;

	// counting reference
	static size_t number_of_instances;
};


#endif /* FFPOPSIM_HIGHD_H_ */
