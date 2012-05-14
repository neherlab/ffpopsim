/**
 * @file popgen_highd.h
 * @brief Header file for high-dimensional simulations
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-19
 *
 * HP_VERBOSE: degree of verbosity of haploid_clone. Levels:
 * - 0: no messages
 * - 1: most messages (enter/exit function)
 * - 2: all messages
 */
#ifndef POPGEN_HIGHD_H_
#define POPGEN_HIGHD_H_
#include "popgen.h"

#define HCF_MEMERR -131545
#define HCF_BADARG -131546
#define HCF_VERBOSE 0

using namespace std;


/**
 * @brief Trait coefficient for a set of loci.
 *
 * This struct is used in the hypercube_function class for saving trait coefficients. See also hypercube_function::add_coefficient.
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
	coeff_single_locus_t(double value_in, int locus_in){
		value=value_in;
		locus=locus_in;
	}
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
class hypercube_function
{
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
	hypercube_function();
	hypercube_function(int dim_in, int s=0);
	virtual ~hypercube_function();
	int set_up(int dim_in,  int s=0);

	// methods
	unsigned int get_dim(){return dim;}
	unsigned int get_seed() {return seed;};
	double get_func(boost::dynamic_bitset<> *genotype);
	double get_additive_coefficient(int locus);
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
	clone_t(int n_traits){trait.resize(n_traits);}
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
class haploid_clone {
public:
	// genotype to traits maps, which in turn are used in the trait-to-fitness map
	hypercube_function *trait;

	// pointers to the current population and the temporary one
	vector <clone_t> *current_pop;
	vector <clone_t> *new_pop;

	// construction / destruction
	haploid_clone();
	virtual ~haploid_clone();
	virtual int set_up(int N_in, int L_in,  int rng_seed=0, int number_of_traits=1);

	// initialization
	int init_genotypes(int n_o_genotypes=-1);
	int init_frequencies(double *nu, int n_o_genotypes=0);

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

	// modify population
	void add_genotypes(boost::dynamic_bitset<> newgt, int n);
	int add_fitness_coefficient(double value, vector <int> loci, int traitnumber=0){return trait[traitnumber].add_coefficient(value, loci);}
	void clear_fitness_function(){for(int t=0; t<number_of_traits; t++){trait[t].coefficients_single_locus.clear(); trait[t].coefficients_epistasis.clear();}}
	void flip_single_locus(unsigned int clonenum, int locus);

	// evolution
	int evolve(int gen=1);	
	int bottleneck(int size_of_bottleneck);

	// update traits and fitness and calculate statistics
	void calc_stat();

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
	vector <int> random_sample;
	void produce_random_sample(int size=1000);

	// population parameters
	int number_of_traits;
	int number_of_individuals_max;		//maximal population size in terms of memory allocated to hold genotypes
	int population_size;				//actual population size
	int scratch;				//variable by how much the memory for offsprings exceeds the number_of_individuals (1+scratch)*..
	int generation;
	int number_of_loci;			//total number of loci

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
};



#endif /* POPGEN_HIGHD_H_ */
