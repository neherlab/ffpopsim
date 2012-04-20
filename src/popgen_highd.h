/**
 * @file popgen_highd.h
 * @brief Header file for high-dimensional simulations
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-19
 */
#ifndef POPGEN_HIGHD_H_
#define POPGEN_HIGHD_H_

#define HCF_MEMERR -131545
#define HCF_BADARG -131546
#define HCF_VERBOSE 0

using namespace std;


/**
 * @brief Trait coefficient for a set of loci.
 *
 * This class is used in the hypercube_function class for saving trait coefficients. See also hypercube_function::add_coefficient.
 *
 * Note: words and bits are deprecated.
 */
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

/**
 * @brief Trait coefficient for a single locus.
 *
 * Note: word and bit are deprecated. For this reason, this class is actually equivalent
 * to the index_value_pair struct.
 */
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

public:
	int dim;
	double hypercube_mean;
	vector <coeff_single_locus> coefficients_single_locus;
	vector <coeff> coefficients_epistasis;
	double epistatic_std;
	int rng_offset;

	// memory management
	bool hcube_allocated;
	bool mem;
	int allocate_mem();
	int free_mem();

	// setting up
	hypercube_function();
	hypercube_function(int dim_in, int s=0);
	~hypercube_function();
	int set_up(int dim_in,  int s=0);

	// methods
	unsigned int get_dim(){return dim;}
	unsigned int get_seed() {return seed;};
	double get_func(boost::dynamic_bitset<> *gt);
	double get_additive_coefficient(int locus);
	int set_additive_coefficient(double value, int locus, int expected_locus=-1);
	int add_coefficient(double value, vector <int> loci);
	int set_random_epistasis_strength(double sigma);
};



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
private:
	int number_of_traits;
	int number_of_individuals;		//maximal population size in terms of memory allocated to hold genotypes
	int pop_size;				//actual population size
	int target_pop_size;			//target (average) population size
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

	boost::dynamic_bitset<> reassortment_pattern();
	boost::dynamic_bitset<> crossover_pattern();

	int flip_single_locus(int locus);
	void flip_single_locus(int individual, int locus);
	void shuffle_genotypes();
	int swap_populations();
	int add_recombinants();
	int new_generation();
	void update_fitness();

protected:
	gsl_rng* evo_generator;
	gsl_rng* label_generator;
	int seed;

public:
	hypercube_function *trait;		// genotype to fitness map
	vector <gt> *current_pop;		// genotypes is a [number_of_individuals*number_of_words] array storing the alleles of each individual
	vector <gt> *new_pop;			// newgenotypes is a [number_of_individuals*number_of_words] array storing the alleles of each individual
	vector <gt> pop2;			// genotypes is a [number_of_individuals*number_of_words] array storing the alleles of each individual
	vector <gt> pop1;			// newgenotypes is a [number_of_individuals*number_of_words] array storing the alleles of each individual

	// constructors/destructors
	haploid_clone();
	haploid_clone(int N_in,int L,  int rng_seed=0, int number_of_traits=1);
	virtual ~haploid_clone();
	int set_up(int N_in,int L,  int rng_seed=0, int number_of_traits=1);

	// initialization
	void set_target_pop_size(int tgps){target_pop_size=tgps;}
	void set_mutation_rate(double mu){mutation_rate=mu;}
	void set_outcrossing_probability(double r){outcrossing_probability=r;}
	void set_fixed_rec_rate(double r){outcrossing_probability=r;}		// deprecated
	void set_crossover_rate(double c){crossover_rate=c;}
	void set_recombination_model(int c) {recombination_model=c;}
	int init_genotypes(int n_o_genotypes=-1);
	int init_genotypes(double *nu, int n_o_genotypes=0);
	void set_circular(bool c) {circular=c;}
	void produce_random_sample(int size);

	// modify population
	void add_genotypes(boost::dynamic_bitset<> newgt,  int n);
	int add_fitness_coefficient(double value, vector <int> loci, int traitnumber=0){return trait[traitnumber].add_coefficient(value, loci);}
	void clear_fitness_function(){for (int t=0; t<number_of_traits; t++){trait[t].coefficients_single_locus.clear(); trait[t].coefficients_epistasis.clear();}}

	// evolve
	int evolve();
	int evolve(int g);
	int bottleneck(int size_of_bottleneck);
	void mutate();
	int select_gametes();
	int recombine(int parent1, int parent2);
	int recombine_crossover(int parent1, int parent2, int ng);
	double chemical_potential();

	// population parameters
	int L(){return number_of_loci;}
	int N() {return pop_size;}
	double mu() {return mutation_rate;}
	int get_pop_size() {return pop_size;}
	int get_target_pop_size() {return target_pop_size;}
	int get_generation(){return generation;}

	// readout
	int random_clone(int size=1000);
	string get_genotype_string(int i);
	void calc_stat();
	void calc_fit();
	void calc_fitness_stat();
	void calc_trait_stat();
	void calc_individual_fitness(gt *tempgt);
	void calc_everybodies_traits();
	virtual void calc_fitness_from_traits(gt *tempgt){tempgt->fitness = tempgt->trait[0];}	// this must be virtual, because the fitness landscape on the (genotype x phenotype) space can be wild
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
	double get_fitness(int n) {return (*current_pop)[n].fitness;}
	double get_trait(int n, int t) {return (*current_pop)[n].trait[t];}
	double get_max_fitness();
	int print_allele_frequencies(ostream &out);
	int read_ms_sample(istream &gts, int skip_locus, int multiplicity);
	int read_ms_sample_sparse(istream &gts, int skip_locus, int multiplicity, int distance);
};


#endif /* POPGEN_HIGHD_H_ */
