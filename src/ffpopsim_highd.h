/**
 * @file popgen_highd.h
 * @brief Header file for high-dimensional simulations
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-04-19
 *
 * Copyright (c) 2012-2013, Richard Neher, Fabio Zanini
 * All rights reserved.
 *
 * This file is part of FFPopSim.
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
 * This struct is used in the hypercube_highd class for saving trait coefficients.
 * See also hypercube_highd::add_coefficient.
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
	// random number generator
	gsl_rng *rng;
	unsigned int seed;

	// memory management
	bool hcube_allocated;
	bool mem;
	int allocate_mem();
	int free_mem();

	// iterators
	vector<coeff_single_locus_t>::iterator coefficients_single_locus_iter;
	vector<coeff_t>::iterator coefficients_epistasis_iter;

	// static array of single locus coefficients (for performance reasons)
	vector<double> coefficients_single_locus_static;

public:
        // random number generator
	int rng_offset;

        // attributes
	int dim;
	double hypercube_mean;
	double epistatic_std;
	vector <coeff_single_locus_t> coefficients_single_locus;
	vector <coeff_t> coefficients_epistasis;

	// setting up
	hypercube_highd();
	hypercube_highd(int dim_in, int s=0);
	virtual ~hypercube_highd();
	int set_up(int dim_in,  int s=0);

	// get methods
	unsigned int get_dim(){return dim;}
	unsigned int get_seed() {return seed;};
	double get_func(boost::dynamic_bitset<>& genotype);
	double get_additive_coefficient(int locus);
	double get_func_diff(boost::dynamic_bitset<>& genotype1, boost::dynamic_bitset<>& genotype2, vector<int> &diffpos);

	// change the hypercube
	void reset();
	void reset_additive();
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
#define HP_VERY_NEGATIVE -1e15

// Error Codes
#define HP_BADARG -879564
#define HP_MEMERR -986465
#define HP_EXPLOSIONWARN 4
#define HP_EXTINCTERR 5
#define HP_NOBINSERR 6
#define HP_WRONGBINSERR 7
#define HP_RUNTIMEERR 8

/**
 * @brief clone with a single genotype and a vector of phenotypic traits.
 *
 * Note: it uses dynamic bitsets because they require little memory.
 */
struct clone_t {
	boost::dynamic_bitset<> genotype;
	vector<double> trait;
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


/*
 *	@brief a class that implements a rooted tree to store genealogies
 *
 *	Nodes and edges are stored as maps with a key that holds the age (rather the time) when the node lived
 *	and the index in the population at that time. The nodes themselves are sufficient to reconstruct the tree
 *	since they contain keys of parents and children
 */

#ifndef rooted_tree_H_
#define rooted_tree_H_
#define RT_VERBOSE 0
#define RT_VERYLARGE 10000000
#define RT_CHILDNOTFOUND -35343
#define RT_NODENOTFOUND -35765
#define RT_LOCUSNOTFOUND -35762
#define RT_FITNESS_MISSING -35722
#define RT_CROSSOVER_MISSING -35721
#define RT_SEGMENT_MISSING -35720
#define RT_ERROR_PARSING 1

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <gsl/gsl_histogram.h>

using namespace std;

struct tree_key_t {
	int index;
	int age;
	bool operator==(const tree_key_t &other)  {return (age == other.age) && (index == other.index);}
	bool operator!=(const tree_key_t &other)  {return (age != other.age) || (index != other.index);}
	bool operator<(const tree_key_t &other) const {
                if(age < other.age) return true;
                else if (age > other.age) return false;
                else { return (index<other.index); }
        }
	bool operator>(const tree_key_t &other) const {
                if(age > other.age) return true;
                else if (age < other.age) return false;
                else { return (index>other.index); }
        }
        tree_key_t(int index=0, int age=0) : index(index), age(age) {};
};

struct step_t {
	int pos;
	int step;
	bool operator<(const step_t &other) const {
                if(pos < other.pos) return true;
                else return false;
        }
	bool operator>(const step_t &other) const {
                if(pos > other.pos) return true;
                else return false;
        }
	bool operator==(const step_t &other) const {
                if(pos == other.pos) return true;
                else return false;
        }
        step_t(int pos=0, int step=0) : pos(pos), step(step) {};
};

struct node_t {
	tree_key_t parent_node;
	tree_key_t own_key;
	list < tree_key_t > child_edges;
	double fitness;
	vector <step_t> weight_distribution;
	int number_of_offspring;
	int clone_size;
	int crossover[2];
};

struct edge_t {
	tree_key_t parent_node;
	tree_key_t own_key;
	int segment[2];
	int length;
	int number_of_offspring;
};

struct poly_t {
	int birth;
	int sweep_time;
	double effect;
	double fitness;
	double fitness_variance;
	poly_t(int b=0, int age=0, double e=0, double f=0, double fvar=0) :
                birth(b), sweep_time(age), effect(e), fitness(f), fitness_variance(fvar) {};
};


class rooted_tree {
public:
	map < tree_key_t , edge_t > edges;
	map < tree_key_t , node_t > nodes;
	vector <tree_key_t> leafs;
	tree_key_t root;
	tree_key_t MRCA;

	rooted_tree();
	virtual ~rooted_tree();
	void reset();
	void add_generation(vector <node_t> &new_generation, double mean_fitness);
	int add_terminal_node(node_t &newNode);
	tree_key_t erase_edge_node(tree_key_t to_be_erased);
	tree_key_t bridge_edge_node(tree_key_t to_be_bridged);
	int external_branch_length();
	int total_branch_length();
	int ancestors_at_age(int age, tree_key_t subtree_root, vector <tree_key_t> &ancestors);
	int update_leaf_to_root(tree_key_t leaf);
	void update_tree();
	int calc_weight_distribution(tree_key_t subtree_root);
	void SFS(gsl_histogram *sfs);
	tree_key_t get_MRCA(){return MRCA;};
	int erase_child(map <tree_key_t,node_t>::iterator Pnode, tree_key_t to_be_erased);
	int delete_extra_children(tree_key_t subtree_root);
	int delete_one_child_nodes(tree_key_t subtree_root);
	bool check_node(tree_key_t node);
	int check_tree_integrity();
	void clear_tree();

        // print tree or subtrees
	string print_newick();
	string subtree_newick(tree_key_t root);
	string print_weight_distribution(tree_key_t node_key);
	int read_newick(string newick_string);

        // construct subtrees
	int construct_subtree(vector <tree_key_t> subtree_leafs, rooted_tree &other);

private:
	static int parse_label(std::string label, int *index, int *clone_size, int *branch_length);
	int parse_subtree(tree_key_t &parent_key, std::string &tree_s);

};

#endif /* rooted_tree_H_ */

/*
 * @brief short wrapper class that handles trees at different places in the genome
 *
 * the class contains a vector of rooted_tree instances that hold the genealogy
 *  in different places. In addition, there is a rooted_tree called subtree
 *  that is used on demand
 *
 *  Created on: Oct 14, 2012
 *      Author: richard
 */

#ifndef MULTILOCUSGENEALOGY_H_
#define MULTILOCUSGENEALOGY_H_
class multi_locus_genealogy {
public:
	vector <int> loci;				//vector of loci (positions on a genome) whose genealogy is to be tracked
	vector <rooted_tree> trees;                     //vector of rooted trees (one per locus)
	vector < vector < node_t > > newGenerations;	//used by the evolving class to store the new generation

	multi_locus_genealogy();
	virtual ~multi_locus_genealogy();
	void track_locus(int new_locus);
	void reset(){loci.clear(); trees.clear();newGenerations.clear();}
	void reset_but_loci(){for(unsigned int i=0; i<loci.size(); i++){trees[i].reset();newGenerations[i].clear();}}
	void add_generation(double baseline);
	int extend_storage(int n);
};
#endif /* MULTILOCUSGENEALOGY_H_ */


/**
 * @brief Population class for high-dimensional simulations.
 *
 * This class is the main object storing the state of and enabling the manipulation of populations with long genomes (\f$L\f$) larger than 20.
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

	// construction / destruction
	haploid_highd(int L=0, int rng_seed=0, int number_of_traits=1, bool all_polymorphic=false);
	virtual ~haploid_highd();

        // the population
	vector <clone_t> population;

	// population parameters (read/write)
	int carrying_capacity;			// carrying capacity of the environment (pop size)
	double outcrossing_rate;		// probability of having sex
	double crossover_rate;			// rate of crossover during sex
	int recombination_model;		//model of recombination to be used
	bool circular;				//topology of the chromosome
	double growth_rate;			//growth rate for bottlenecks and the like

        // mutation rate (only if not all_polymorphic)
        double get_mutation_rate(){return mutation_rate;}
        void set_mutation_rate(double m){
        if(all_polymorphic){
                if(HP_VERBOSE) cerr<<"Cannot set the mutation rate with all_polymorphic."<<endl;
                throw HP_BADARG;
        } else mutation_rate=m;}

        // pseudo-infinite site model
        bool is_all_polymorphic(){return all_polymorphic;}
        vector<poly_t> get_polymorphisms(){return polymorphism;}
        vector<poly_t> get_fixed_mutations(){return fixed_mutations;}
        vector<int> get_number_of_mutations(){return number_of_mutations;}

	// population parameters (read only)
	int L(){return number_of_loci;}
	int get_number_of_loci(){return number_of_loci;}
	int N(){return population_size;}
	int get_population_size() {return population_size;}
	int get_generation(){return generation;}
	void set_generation(int g){generation = g;}
	int get_number_of_clones(){return number_of_clones;}
	int get_number_of_traits(){return number_of_traits;}
	double get_participation_ratio(){return participation_ratio;}

	// initialization
	int set_allele_frequencies(double* frequencies, unsigned long N);
	int set_genotypes_and_ancestral_state(vector <genotype_value_pair_t> gt, vector <int> anc_state);
	int set_genotypes(vector <genotype_value_pair_t> gt);
	int set_wildtype(unsigned long N);
	int track_locus_genealogy(vector <int> loci);

	// modify population
	void add_genotype(boost::dynamic_bitset<> genotype, int n=1);

	// modify traits
	int add_trait_coefficient(double value, vector <int> loci, int t=0){return trait[t].add_coefficient(value, loci);}
	void clear_trait(int t=0){if(t >= number_of_traits) throw (int)HP_BADARG; else trait[t].reset();}
	void clear_traits(){for(int t=0; t<number_of_traits; t++){trait[t].reset();}}
	void set_random_trait_epistasis(double epistasis_std,int traitnumber=0){trait[traitnumber].epistatic_std=epistasis_std;}

	// modify fitness (shortcuts: they only make sense if number_of_traits=1)
	int add_fitness_coefficient(double value, vector <int> loci){if(number_of_traits>1) throw (int)HP_BADARG; return add_trait_coefficient(value, loci, 0);}
	void clear_fitness(){if(number_of_traits>1){if(HP_VERBOSE) cerr<<"What do you mean by fitness?"<<endl; throw (int)HP_BADARG;} clear_traits();}
	void set_random_epistasis(double epistasis_std){if(number_of_traits>1){if(HP_VERBOSE) cerr<<"Please use set_random_trait_epistasis."<<endl; throw (int)HP_BADARG;} trait[0].epistatic_std=epistasis_std;}

	// evolution
	int evolve(int gen=1);	
	int bottleneck(int size_of_bottleneck);
	unsigned int flip_single_locus(int locus);

	// update traits and fitness and calculate statistics
	void calc_stat();
	void unique_clones();
        vector <int> get_nonempty_clones();

	// readout
	// Note: these functions are for the general public and are not expected to be
	// extremely fast. If speed is a major concern, consider subclassing and working
	// with protected methods.

	// random clones
	int random_clone();
	int random_clones(unsigned int n_o_individuals, vector <int> *sample);

	// genotype readout
	string get_genotype_string(unsigned int i){string gts; boost::to_string(population[i].genotype, gts); return gts;}
	int distance_Hamming(unsigned int clone1, unsigned int clone2, vector <unsigned int *> *chunks=NULL, unsigned int every=1){return distance_Hamming(population[clone1].genotype, population[clone2].genotype, chunks, every);}
	int distance_Hamming(boost::dynamic_bitset<> gt1, boost::dynamic_bitset<> gt2, vector<unsigned int *> *chunks=NULL, unsigned int every=1);
	stat_t get_diversity_statistics(unsigned int n_sample=1000);
	stat_t get_divergence_statistics(unsigned int n_sample=1000);

	// allele frequencies
	double get_allele_frequency(int l) {if (!allele_frequencies_up_to_date){calc_allele_freqs();} return allele_frequencies[l];}
	double get_derived_allele_frequency(int l) {if (ancestral_state[l]) {return 1.0-get_allele_frequency(l);} else {return get_allele_frequency(l);}}
	bool get_ancestral_state(int l) {return ancestral_state[l];}

	double get_pair_frequency(int locus1, int locus2);
	vector <double> get_pair_frequencies(vector < vector <int> > *loci);
	double get_chi(int l) {return 2 * get_allele_frequency(l) - 1;}
	double get_derived_chi(int l) {return 2 * get_derived_allele_frequency(l) - 1;}
	double get_chi2(int locus1, int locus2){return get_moment(locus1, locus2)-get_chi(locus1)*get_chi(locus2);}
	double get_LD(int locus1, int locus2){return 0.25 * get_chi2(locus1, locus2);}
	double get_moment(int locus1, int locus2){return 4 * get_pair_frequency(locus1, locus2) + 1 - 2 * (get_allele_frequency(locus1) + get_allele_frequency(locus2));}

	// fitness/phenotype readout
	void set_trait_weights(double *weights){for(int t=0; t<number_of_traits; t++) trait_weights[t] = weights[t];}
	double get_trait_weight(int t){return trait_weights[t];}
	double get_fitness(int n) {calc_individual_fitness(population[n]); return population[n].fitness;}
	int get_clone_size(int n) {return population[n].clone_size;}
	double get_trait(int n, int t=0) {calc_individual_traits(population[n]); return population[n].trait[t];}
	vector<coeff_t> get_trait_epistasis(int t=0){return trait[t].coefficients_epistasis;}
	stat_t get_fitness_statistics() {update_fitness(); calc_fitness_stat(); return fitness_stat;}
	stat_t get_trait_statistics(int t=0) {calc_trait_stat(); return trait_stat[t];}
	double get_trait_covariance(int t1, int t2) {calc_trait_stat(); return trait_covariance[t1][t2];}
	double get_max_fitness() {return fitness_max;}
	void update_traits();
	void update_fitness();

	// histograms
	int get_divergence_histogram(gsl_histogram **hist, unsigned int bins=10, vector <unsigned int *> *chunks=NULL, unsigned int every=1, unsigned int n_sample=1000);
	int get_diversity_histogram(gsl_histogram **hist, unsigned int bins=10, vector <unsigned int *> *chunks=NULL, unsigned int every=1, unsigned int n_sample=1000);
	int get_fitness_histogram(gsl_histogram **hist, unsigned int bins=10, unsigned int n_sample=1000);

	// stream I/O
	int print_allele_frequencies(ostream &out);
	int read_ms_sample(istream &gts, int skip_locus, int multiplicity);
	int read_ms_sample_sparse(istream &gts, int skip_locus, int multiplicity, int distance);

        // genealogy
	multi_locus_genealogy genealogy;

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
	int generation;
	int number_of_clones;
	double mutation_rate;			// rate of mutation per locus per generation

	// evolution
	int mutate();
	int select_gametes();
	double relaxation_value();
	double get_logmean_expfitness();	// Log of the population exp-average of the fitness: log[<exp(F)>_{population}]
	
	unsigned int flip_single_locus(unsigned int clonenum, int locus);
	void shuffle_genotypes();
	int new_generation();

	// clone structure
	double participation_ratio;
	int partition_cumulative(vector <unsigned int> &partition_cum);
	int provide_at_least(int n);
	int last_clone;

	// allele_frequencies
	bool allele_frequencies_up_to_date;
	double *allele_frequencies;
	double *gamete_allele_frequencies;
	double *chi1;				//symmetric allele frequencies
	double **chi2;				//symmetric two locus correlations
	bool all_polymorphic;                   // switch that makes sure every locus is polymorphic in an infinite alleles model when mutation rate is 0
	vector <int> ancestral_state;	//vector, that for each locus keeps track of the ancestral state. by default, all zero
	vector <poly_t> polymorphism;	//vector, that keeps track when an allele was introduced on which background. Only needed in an infinite alleles model
	vector <poly_t> fixed_mutations;	//vector to store all fixed mutations
	vector <int> number_of_mutations;	//vector to store the number of mutations introduced each generation
	void calc_allele_freqs();

	// recombination details
	double outcrossing_rate_effective;
	int *genome;				//Auxiliary array holding the positions along the genome
	int *crossovers;
	void reassortment_pattern();
	void crossover_pattern();
	vector <int> sex_gametes;		//array holding the indices of gametes
	int add_recombinants();
	int recombine(int parent1, int parent2);
	int recombine_crossover(int parent1, int parent2, int ng);

	// fitness and traits
	double fitness_max;
	stat_t fitness_stat;
	stat_t *trait_stat;
	double **trait_covariance;
	void calc_fitness_stat();
	void calc_trait_stat();
	void calc_individual_traits(clone_t &tempgt);
	void calc_individual_fitness(clone_t &tempgt);
	void calc_individual_traits(int clonenum){calc_individual_traits(population[clonenum]);}
	void calc_individual_fitness(int clonenum){calc_individual_fitness(population[clonenum]);}
	void check_individual_maximal_fitness(clone_t &tempgt){fitness_max = fmax(fitness_max, tempgt.fitness);}
	double get_trait_difference(clone_t &tempgt1, clone_t &tempgt2, vector<int>& diffpos, int traitnum);

	// phenotype-fitness map. By default, a linear map with equal weights is set, but weights can be reset
	double *trait_weights;
	virtual void calc_individual_fitness_from_traits(clone_t &tempgt);
	virtual void calc_individual_fitness_from_traits(int clonenum) {calc_individual_fitness_from_traits(population[clonenum]);}
	void add_clone_to_genealogy(int locus, int dest, int parent, int left, int right, int cs, int n);
	bool track_genealogy;

private:
	// Memory management is private, subclasses must take care only of their own memory
	bool mem;
	bool cumulants_mem;
	int allocate_mem();
	int free_mem();

	// These two vectors are used to recycle dead clones
	vector <int> available_clones;
	vector <int> clones_needed_for_recombination;

	boost::dynamic_bitset<> rec_pattern;

	// counting reference
	static size_t number_of_instances;
};

#endif /* FFPOPSIM_HIGHD_H_ */
