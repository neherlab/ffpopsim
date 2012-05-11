// vim: tabstop=8:softtabstop=8:shiftwidth=8:noexpandtab
/*
 * haploid_clone.cpp
 *
 *  Created on: Aug 28, 2008
 *      Author: neher
 *    Modified: Dec 13, 2011
 *   Committer: Fabio Zanini
 */
#include <math.h>
#include "popgen_highd.h"

/**
 * @brief Default constructor.
 *
 * The population objects are initialized, but neither the random generator is initialized nor the memory
 * for traits, gametes etc. is allocated. Please call set_up() before using the class.
 *
 * Note: The sequence is assumed to be linear (not circular). You can change this by hand if you wish so.
 */
haploid_clone::haploid_clone() {
	current_pop = &current_pop_vector;
	new_pop = &new_pop_vector;
	mem=false;
	cumulants_mem=false;
	circular=false;
}

/**
 * @brief Destructor.
 *
 * Memory is released here.
 */
haploid_clone::~haploid_clone() {
	if (mem) free_mem();
}

/**
 * @brief Construct a population with certain parameters.
 *
 * @param N_in number of individuals
 * @param L_in length of the genome
 * @param rng_seed seed for the random number generator. If this is 0, time(NULL)+getpid() is used.
 * @param n_o_traits number of phenotypic traits (including fitness). Must be \f$\geq1\f$.
 *
 * @returns zero if successful, error codes otherwise
 *
 * Note: memory allocation is also performed here, via the allocate_mem function.
 */
int haploid_clone::set_up(int N_in, int L_in,  int rng_seed, int n_o_traits)
{
	boost::dynamic_bitset<> temp;
	if (N_in<1 or L_in <1 or n_o_traits<1)
	{
		cerr <<"haploid_clone::set_up(): Bad Arguments! All of N, L and the number of traits must be larger or equal one.\n";
		return HP_BADARG;
	}

	if (8*sizeof(long int)!=temp.bits_per_block) {
		cerr <<"haploid_clone requires sizeof(long int) to be equal to bits_per_block of boost::dynamic_bitset";
		return HP_BADARG;
	}

	number_of_traits=n_o_traits;
	number_of_loci=L_in;
	number_of_individuals_max=2*N_in;
	carrying_capacity=N_in;

	//In case no seed is provided use current second and add process ID
	if (rng_seed==0){
		seed=time(NULL)+getpid();
	}else{seed=rng_seed;}

	mem=false;
	int err=allocate_mem();
	return err;
}

/**
 * @brief Allocate all the necessary memory, initialze the RNG.
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_clone::allocate_mem()
{
	if (mem)
	{
		cerr <<"haploid_clone::allocate_mem(): memory already allocated, freeing and reallocating ...!\n";
		free_mem();
	}

	if (HP_VERBOSE) cerr <<"haploid_clone::allocate_mem(): genotypes: "<<number_of_individuals_max<<"  loci: "<<number_of_loci<<endl;

	//Random number generator
	evo_generator=gsl_rng_alloc(RNG);
	gsl_rng_set(evo_generator, seed);
	//allocate all the memory
	genome = new int [number_of_loci+1];					// aux array holding range(0,number_of_loci) used to draw crossover points
	for (int i=0; i<number_of_loci; i++) genome[i]=i;
	crossovers= new int [number_of_loci];					// aux array holding crossover points

	if (HP_VERBOSE) cerr <<"allele frequencies...";
	allele_frequencies =new double [number_of_loci];
	gamete_allele_frequencies =new double [number_of_loci];			//allele frequencies after selection

	trait = new hypercube_function [number_of_traits];			//genotype trait function
	trait_stat = new stat_t [number_of_traits];				//structure holding trait statistics
	trait_covariance = new double* [number_of_traits];
	//initialize trait functions
	for (int t=0; t<number_of_traits;t++){
		trait[t].set_up(number_of_loci, gsl_rng_uniform_int(evo_generator, 1<<20));
		trait_covariance[t]=new double [number_of_traits];
	}

	mem=true;								//set memory flag to true
	if (HP_VERBOSE) cerr <<"done.\n";
	generation=0;
	return 0;
}

/**
 * @brief Releases memory during ckass destruction.
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_clone::free_mem()
{
	if (!mem)
	{
		cerr <<"haploid_clone::free_mem(): No memory allocated!\n";
		return HP_BADARG;
	}
	delete [] allele_frequencies;
	delete [] gamete_allele_frequencies;
	delete [] trait;
	mem=false;
	return 0;
}

/**
 * @brief Initialize population with a certain number of individuals and fixed allele frequencies.
 *
 * @param nu target allele frequencies
 * @param n_o_genotypes number of individuals to be created
 *
 * @returns zero if successful, error codes otherwise
 *
 * Note: when this function is used to initialize the population, it is likely that the fitness distribution
 * has a very large width. In turn, this can result in an immediate and dramatic drop in diversity within the
 * first few generations. Please check fitness statistics before starting the evolution if this worries you.
 */
int haploid_clone::init_genotypes(double* nu, int n_o_genotypes)
{
	int i, locus;
	if (!mem)
	{
		cerr <<"haploid_clone::init_genotypes(): Allocate and set up first!\n";
		return HP_MEMERR;
	}
	if (n_o_genotypes<0 or n_o_genotypes>=number_of_individuals_max)
	{
		cerr <<"haploid_clone::init_genotypes(): number of genotypes has to be positive and smaller than max! got: "<<n_o_genotypes<<"!\n";
		return HP_BADARG;
	}
	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes(double* nu, int n_o_genotypes) ....";
	generation=0;
	int ngt=0;
	// if number of genotypes to be drawn is not specified, use default carrying capacity
	if (n_o_genotypes==0) ngt=carrying_capacity;
	else ngt=n_o_genotypes;

	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes(double* nu, int n_o_genotypes) reset current\n";
	current_pop->clear();	//reset the current population
	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes(double* nu, int n_o_genotypes) reset random sample\n";
	random_sample.clear();	//and the random sample
	boost::dynamic_bitset<> tempgt(number_of_loci);
	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes(double* nu, int n_o_genotypes) add "<<ngt<<" genotypes of length "<<number_of_loci<<"\n";
	for (i=0; i<ngt; i++) {
		tempgt.reset();
		for(locus=0; locus<number_of_loci; locus++){	//set all loci
			if (gsl_rng_uniform(evo_generator)<nu[locus])	tempgt.set(locus);
		}
		add_genotypes(tempgt,1);	//add genotype with multiplicity 1
		if (HP_VERBOSE) cerr <<i<<endl;
	}

	generation=0;
	calc_stat();

	if (HP_VERBOSE) cerr <<"done.\n";
	return 0;
}



/**
 * @brief Initialize with a single genotypic clone (00...0).
 *
 * @param n_o_genotypes size of the clone. If not chosen, use carrying_capacity.
 *
 * @returns 0 if successful, nonzero otherwise.
 *
 * This is the typical initialization function if evolution from a single ancestor is modeled.
 */
int haploid_clone::init_genotypes(int n_o_genotypes)
{
	trait[0].coefficients_epistasis.clear();

	if (!mem)
	{
		cerr <<"haploid_clone::init_genotypes(): Allocate and set up first!\n";
		return HP_MEMERR;
	}
	if (n_o_genotypes>number_of_individuals_max)
	{
		cerr <<"haploid_clone::init_genotypes(): number of genotypes has to be smaller than max! got: "<<n_o_genotypes<<"!\n";
		return HP_BADARG;
	}
	population_size=0;
	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes(int n_o_genotypes) with population size "<<n_o_genotypes<<"...";
	current_pop->clear();
	random_sample.clear();
	boost::dynamic_bitset<> tempgt(number_of_loci);
	tempgt.reset();					//set all bits to zero
	if (n_o_genotypes>=0){				//add n_o_genotypes copies of this genotype
		add_genotypes(tempgt,n_o_genotypes);
	}else{
		add_genotypes(tempgt,carrying_capacity);
	}

	generation=0;
	calc_stat();

	if (HP_VERBOSE) cerr <<"done."<<endl;
	return 0;
}


/**
 * @brief calculate and store allele frequencies.
 *
 * Note: the allele frequencies are available in the allele_frequencies attribute.
 */
void haploid_clone::calc_allele_freqs()
{
	if (HP_VERBOSE) cerr<<"haploid_clone::calc_allele_freqs()...";
	int locus,cs;
	population_size=0;
	for (locus=0; locus<number_of_loci; allele_frequencies[locus++] = 0);
	//loop over all clones
	for (unsigned int i=0; i<current_pop->size(); i++)
	{
		cs=(*current_pop)[i].clone_size;
		for (locus=0; locus<number_of_loci; locus++)	//add clone size to allele frequency of clone carries allele
			if ((*current_pop)[i].genotype[locus])
				allele_frequencies[locus] += cs;
		population_size += cs;
	}
	//convert counts into frequencies
	for (locus=0; locus<number_of_loci; allele_frequencies[locus++] /= population_size);
	if (HP_VERBOSE) cerr<<"done.\n";
}

/**
 * @brief Get the joint frequency of two alleles.
 *
 * @param locus1 position of the first allele.
 * @param locus2 position of the second allele.
 *
 * @returns the joint frequency of the two alleles.
 */
double haploid_clone::get_pair_frequency(int locus1, int locus2)
{
	if (HP_VERBOSE) cerr<<"haploid_clone::get_pair_frequency()...";
	double frequency = 0;
	for (unsigned i=0; i<current_pop->size(); i++)
		if ((*current_pop)[i].genotype[locus1] and (*current_pop)[i].genotype[locus2])
		        frequency += (*current_pop)[i].clone_size;
	frequency /= population_size;
	if (HP_VERBOSE) cerr<<"done.\n";
	return frequency;
}

/**
 * @brief Get the joint frequency of two alleles, for a vector of allele pairs.
 *
 * @param loci pointer to a vector of allele pairs. Each element of loci must be a vector of length 2.
 *
 * @returns vector of joint frequencies.
 */
vector <double>  haploid_clone::get_pair_frequencies(vector < vector <int> > *loci)
{
	if (HP_VERBOSE) cerr<<"haploid_clone::get_pair_frequencies()...";
	unsigned int pair;

	vector <double> freq (loci->size(), 0.0);
	for (pair = 0; pair < loci->size(); pair++)
		freq[pair] = get_pair_frequency((*loci)[pair][0], (*loci)[pair][1]);
	if (HP_VERBOSE) cerr<<"done.\n";
	return freq;
}


/**
 * @brief Evolve for some generations under the specified conditions.
 *
 * The order of steps performed is
 * 1. mutation
 * 2. selection
 * 3. recombination
 * but should not be important except for extremely high selection coefficients and/or very short times,
 * for which the discrete nature of the Fisher-Wright model becomes relevant.
 *
 * @param gen number of generations.
 *
 * @returns sum of error codes of the single evolution steps. It is therefore a multiple of gen.
 *
 * Note: if an error in encountered, evolution is stopped after the function that created the problem.
 * Typical errors include extinction or the opposite, offspring explosion.
 */
int haploid_clone::evolve(int gen){
	if (HP_VERBOSE) cerr<<"haploid_clone::evolve(int gen)...";
	int err=0, gtemp = 0, btn = 0;

	// calculate an effective outcrossing rate to include the case of very rare crossover rates.
	// Since a recombination without crossovers is a waste of time, we scale down outcrossing probability
	// and scale up crossover rate so that at least one crossover is guaranteed to happen.
	if (recombination_model==CROSSOVERS)
		outcrossing_rate_effective = outcrossing_rate * (1 - exp(- number_of_loci * crossover_rate));
	else
		outcrossing_rate_effective = outcrossing_rate;

	// evolve cycle
	while((err == 0) && (gtemp < gen)) {
		if (HP_VERBOSE) cerr<<"generation "<<(generation+gtemp)<<endl;
		random_sample.clear();			//discard the old random sample
		if(err==0) err=mutate();		//mutation step
		if(err==0) err=select_gametes();	//select a new set of gametes (partitioned into sex and asex)
		if(err==0) err=add_recombinants();	//do the recombination between pairs of sex gametes
		if(err==0) err=swap_populations();	//make the new population the current population
		gtemp ++;
	}
	generation+=gtemp;
	if (HP_VERBOSE) {
		if(err==0) cerr<<"done."<<endl;
		else cerr<<"error "<<err<<"."<<endl;
	}
	return err;
}

/**
 * @brief Generate offpring according to fitness (selection) and segregate some for sexual mating. 
 *
 * @returns zero if successful, error codes otherwise
 *
 * Random Poisson offspring numbers are drawn from all parents, proportionally to their fitness.
 *
 * Note: this function needs a few quirks to work properly. The problem has to do with very fit individuals,
 * which would theoretically leave behind a huge number of offspring that immediately fill up the ecological
 * niche. If we follow the rule blindly, we spend a lot of time generating random numbers for any kind of
 * process (mutation, recombination, etc.). Real organisms are limited in their offspring number by physical
 * constraints, which are abstracted away at the level of this library. Thus, we have to use a practical
 * approach. We apply two corrections:
 * 1. the max number of offspring of each clone is limited by MAX_DELTAFITNESS
 * 2. the max population size is reduced below MAX_POPSIZE before exiting the function.
 */
int haploid_clone::select_gametes() {
	if (HP_VERBOSE) cerr<<"haploid_clone::select_gametes()...";
	//determine the current mean fitness, which includes a term to keep the population size constant
	double relaxation = relaxation_value();

	//draw gametes according to parental fitness
	double delta_fitness;
	int os,o, nrec;
	int err=0;
	population_size=0;
	
	//to speed things up, reserve the expected amount of memory for sex gametes and the new population (+10%)
	sex_gametes.reserve(population_size*outcrossing_rate_effective*1.1);
	sex_gametes.clear();
	new_pop->reserve(current_pop->size()*1.1);
	new_pop->clear();
	for (unsigned int i=0; i<current_pop->size(); i++) {
		//poisson distributed random numbers -- mean exp(f)/bar{exp(f)}) (since death rate is one, the growth rate is (f-bar{f})
		if ((*current_pop)[i].clone_size>0){
			//the number of asex offspring of clone[i] is poisson distributed around e^F / <e^F> * (1-r)
			delta_fitness = (*current_pop)[i].fitness - relaxation;
//			if (HP_VERBOSE >= 2) cerr<<i<<": relative fitness = "<<delta_fitness<<", Poisson intensity = "<<((*current_pop)[i].clone_size*exp(delta_fitness)*(1-outcrossing_rate_effective))<<endl;
			os=gsl_ran_poisson(evo_generator, (*current_pop)[i].clone_size*exp(delta_fitness)*(1-outcrossing_rate_effective));

			if (os>0){
				// clone[i] to new_pop with os as clone size
				new_pop->push_back((*current_pop)[i]);
				new_pop->back().clone_size=os;
				population_size+=os;
			}
			//draw the number of sexual offspring, add them to the list of sex_gametes one by one
			if (outcrossing_rate_effective>0){
				nrec=gsl_ran_poisson(evo_generator, (*current_pop)[i].clone_size*exp(delta_fitness)*outcrossing_rate_effective);
				for (o=0; o<nrec; o++) sex_gametes.push_back(i);
			}
		}
	}
	
	if(population_size < 1) {
		err = HP_EXTINCTERR;
		if (HP_VERBOSE) cerr<<"error "<<err<<". The population went extinct!"<<endl;
	}
	if ((HP_VERBOSE) && (err==0)) cerr<<"done."<<endl;
	return err;
}

/**
 * @brief Cause a bottleneck in the population size.
 *
 * @param size_of_bottleneck number of individuals to leave alive (approximate)
 *
 * @returns zero if successful, error codes otherwise
 *
 * The bottleneck is performed as follows. Each clone is traded in for a smaller one with size given by
 * a Poisson random number around the new population size times the clone frequency. This function should
 * therefore almost conserve genotype frequencies, except for rare genotypes that are lost.
 *
 * TODO: this function should accept a gsl random distribution as optional argument for choosing how sharp
 * the bottleneck should be (i.e., how large fluctuations around the expected frequency may be). However,
 * this requires function pointers or templates or lambda functions, and might be a nightmare to code.
 */
int haploid_clone::bottleneck(int size_of_bottleneck) {
	double ostmp;
	unsigned int os;
	int err=0;
	unsigned int old_size = population_size;
	if (HP_VERBOSE) cerr<<"haploid_clone::bottleneck()...";

	population_size = 0;
	new_pop->reserve(size_of_bottleneck*1.1);
	new_pop->clear();

	for(size_t i=0; i < current_pop->size(); i++) {
		ostmp = (*current_pop)[i].clone_size * size_of_bottleneck / double(old_size);
		os=gsl_ran_poisson(evo_generator, ostmp);
		if(os > 0) {
			new_pop->push_back((*current_pop)[i]);
			new_pop->back().clone_size=os;
			population_size+=os;
		}
	}

	if(population_size < 1) err = HP_EXTINCTERR;
	else swap_populations();
	if (HP_VERBOSE) {
		if(err==0) cerr<<"done."<<endl;
		else if(err == HP_EXTINCTERR) cerr<<" The population went extinct!"<<endl;
	}

	return err;
}

/**
 * @brief Mutate random clones at random loci (all loci).
 *
 * Note: for efficiency reasons, if the mutation rate is zero, an if cycle skips the body altogether.
 * The user can therefore call this function safely even if in non-mutating populations, without
 * loss of performance.
 */
int haploid_clone::mutate() {
	if (HP_VERBOSE)	cerr <<"haploid_clone::mutate() ..."<<endl;
	int i, actual_n_o_mutations,locus=0;
	if(mutation_rate) {
		produce_random_sample(number_of_loci*mutation_rate*population_size*2);
		for (locus = 0; locus<number_of_loci; locus++) {
			//draw the number of mutation that are to happen at this locus
			actual_n_o_mutations = gsl_ran_poisson(evo_generator, mutation_rate*population_size);
			//introduce these mutations one by one
			for (i = 0; i < actual_n_o_mutations; ++i) {
				flip_single_locus(random_clone(), locus);
			}
		}
	}
	else if(HP_VERBOSE)
		cerr<<"the mutation rate is zero...";
	if (HP_VERBOSE)	cerr <<"done."<<endl;;
	return 0;
}

/**
 * @brief Flip a spin at a specific locus in random individual.
 *
 * @param locus position of the locus to flip
 *
 * @returns the (random) flipped individual clone
 *
 * Note: this function calls flip_single_locus(unsigned int clonenum, int locus).
 */
int haploid_clone::flip_single_locus(int locus)
{
	int clonenum = random_clone();
	flip_single_locus(clonenum, locus);
	return clonenum;
}


/**
 * @brief Flip a spin (locus) in individual.,
 *
 * @param clonenum the individual whose locus is being flipped
 * @param locus position of the locus to flip
 *
 * This function creates a new clone and adds it to the population,
 * and assigns it a fitness.
 *
 * Note: in principle, we should look whether an identical clone already exists and,
 * in positive case, add this individual to that clone instead of starting a new one.
 * However, this would take forever, and is thus implemented in a separate function
 * (TODO: which one?).
 */
void haploid_clone::flip_single_locus(unsigned int clonenum, int locus)
{
	//produce new genotype
	clone_t tempgt(number_of_traits);
	tempgt.genotype.resize(number_of_loci);
	tempgt.genotype= (*current_pop)[clonenum].genotype;
	//new clone size == 1, old clone reduced by 1
	tempgt.clone_size=1;
	(*current_pop)[clonenum].clone_size--;
	//flip the locus in new clone and calculate fitness
	tempgt.genotype.flip(locus);
	calc_individual_fitness(&tempgt);
	//add clone to current population
	current_pop->push_back(tempgt);
	if (HP_VERBOSE >= 2) cerr <<"subpop::flip_single_spin(): mutated individual in clone "<<clonenum<<" at locus "<<locus<<endl;
}

/**
 * @brief Pair and mate sexual gametes.
 *
 * @returns zero if successful, nonzero otherwise
 *
 * Using the previously produced list of sex_gametes, pair them at random and mate
 */
int haploid_clone::add_recombinants()
{
	//construct new generation
	int parent1, parent2, n_o_c=1;

	if (sex_gametes.size()>0)
	{
		//sexual offspring -- shuffle the set of gametes to ensure random mating
		gsl_ran_shuffle(evo_generator, &sex_gametes[0], sex_gametes.size(), sizeof(int));
		//if the number of sex_gametes is odd, the last one is unlucky
		for (unsigned int i = 0; i < sex_gametes.size()-1; i+=2) {
			parent1=sex_gametes[i];
			parent2=sex_gametes[i+1];
			//The recombination function stores two new genotypes
			//in new_genotypes[new_population_size] and [new_population_size+1]
			n_o_c=recombine(parent1, parent2);
		}
	}
	return 0;
}

/**
 * @brief Make the the temporary population the current one.
 *
 * After the new population is completely assembled, we have to make it the current
 * population. This is implemented simply by swaping pointers.
 */
int haploid_clone::swap_populations(){
	vector <clone_t> *temp_gt;
	temp_gt=current_pop;
	current_pop=new_pop;
	new_pop=temp_gt;

	return 0;
}


/*
 * recombine two genotypes parent1 and parent2 to produce two new genotypes
 * stored in new_pop at positions ng and ng+1
 *
 */
/**
 * @brief Recombine two genotypes parent1 and parent2 to produce two new genotypes.
 *
 * @param parent1 first parent
 * @param parent2 second parent
 *
 * @returns one (FIXME?)
 *
 * The new genotypes are stored in new_pop at positions ng and ng+1.
 *
 */
int haploid_clone::recombine(int parent1, int parent2) {
	if(HP_VERBOSE >= 2) cerr<<"haploid_clone::recombine(int parent1, int parent2)..."<<endl;

	boost::dynamic_bitset<> rec_pattern;
	//depending on the recombination model, produce a map that determines which offspring
	//inherites which part of the parental genomes
	if (recombination_model==FREE_RECOMBINATION){
		rec_pattern= reassortment_pattern();
	}
	else if (recombination_model==CROSSOVERS)
	{
		rec_pattern= crossover_pattern();
	}
	else {rec_pattern.resize(number_of_loci);}

	//produce two new genoytes
	clone_t offspring1(number_of_traits);
	clone_t offspring2(number_of_traits);
	offspring1.genotype.resize(number_of_loci);
	offspring2.genotype.resize(number_of_loci);

	//assign the genotypes by combining the relevant bits from both parents
	offspring1.genotype=((*current_pop)[parent1].genotype&rec_pattern)|((*current_pop)[parent2].genotype&(~rec_pattern));
	offspring2.genotype=((*current_pop)[parent2].genotype&rec_pattern)|((*current_pop)[parent1].genotype&(~rec_pattern));
	//clone size of new genoytpes is 1 each
	offspring1.clone_size=1;
	offspring2.clone_size=1;

	//Check what's going on
	if(HP_VERBOSE >= 3) {
		cerr<<rec_pattern<<endl;
		cerr<<(*current_pop)[parent1].genotype<<endl;
		cerr<<(*current_pop)[parent2].genotype<<endl;
		cerr<<offspring1.genotype<<endl;
		cerr<<offspring2.genotype<<endl<<endl;
	}

	calc_individual_fitness(&offspring1);
	calc_individual_fitness(&offspring2);

	//add genotypes to new population
	new_pop->push_back(offspring1);
	new_pop->push_back(offspring2);
	population_size+=2;

	if(HP_VERBOSE >= 2) cerr<<"done."<<endl;
	return 1;
}

/**
 * @brief For each clone, recalculate its traits.
 */
void haploid_clone::update_traits() {
	for (unsigned int i=0; i<current_pop->size(); i++) {
		calc_individual_traits(&(*current_pop)[i]);
	}
}

/**
 * @brief For each clone, update fitness assuming traits are already up to date.
 */
void haploid_clone::update_fitness() {
	for (unsigned int i=0; i<current_pop->size(); i++) {
		calc_individual_fitness_from_traits(&(*current_pop)[i]);
	}
}

/**
 * @brief Calculate trait and fitness statistics and allele frequences.
 *
 * Four things are done in a row:
 * 1. traits are updated
 * 2. fitness is updated based on the new traits
 * 3. statistics of traits and fitness are calculated
 * 4. allele freqs are calculated
 *
 * Note: This function is quite expensive. Please use its subblocks
 * whenever possible, and rely on this only when you want to make sure,
 * at the expense of performance, that everything is up to date.
 */
void haploid_clone::calc_stat() {
	update_traits();
	update_fitness();
	calc_trait_stat();
	calc_fitness_stat();
	calc_allele_freqs();
}

/**
 * @brief Calculate traits of the chosen clone.
 *
 * @param tempgt clone whose traits are to be calculated
 */
void haploid_clone::calc_individual_traits(clone_t *tempgt){
	for (int t=0; t<number_of_traits; t++){
		tempgt->trait[t] = trait[t].get_func(&(tempgt->genotype));
	}
}

/**
 * @brief Calculate fitness of a particular clone.
 *
 * @param tempgt clone whose fitness is being calculated
 *
 * Note: this function also updates the traits information for the same clone, because the
 * phenotype is needed to calculate fitness. If you have already calculated the traits,
 * you can rely calc_individual_fitness_from_traits.
 */
void haploid_clone::calc_individual_fitness(clone_t *tempgt) {
	//calculate the new fitness value of the mutant
	calc_individual_traits(tempgt);
	calc_individual_fitness_from_traits(tempgt);
}

/**
 * @brief Choose a certain number of crossover points and produce a 0000111100011111101101
 * crossover pattern
 *
 * @returns crossover pattern
 */
boost::dynamic_bitset<> haploid_clone::crossover_pattern() {
	int n_o_c=0;

	double total_rec = number_of_loci * crossover_rate;
	//TODO this should be poisson conditional on having at least one
	if (total_rec<0.1) n_o_c=1;
	else while (n_o_c==0) n_o_c= gsl_ran_poisson(evo_generator,number_of_loci*crossover_rate);

	//for circular chromosomes make sure there is an even number of crossovers
	if (circular)
	{
		n_o_c*=2;
		n_o_c=(n_o_c<number_of_loci)?n_o_c:number_of_loci;	//make sure there are fewer xovers than loci
		//choose xovers at random from the genome label list
		gsl_ran_choose(evo_generator,crossovers,n_o_c,genome,number_of_loci,sizeof(int));
	}
	else
	{
		n_o_c=(n_o_c<number_of_loci)?n_o_c:(number_of_loci-1);
		gsl_ran_choose(evo_generator,crossovers,n_o_c,(genome+1),number_of_loci-1,sizeof(int));
	}

	bool origin=true;
	boost::dynamic_bitset<> rec_pattern;
	//start with an empty bitset and extend to crossovers[c] with origing =0,1
	for (int c=0;c<n_o_c; c++){
		rec_pattern.resize(crossovers[c], origin);
		origin=!origin; //toggle origin
	}
	//extend to full length
	rec_pattern.resize(number_of_loci, origin);

	return rec_pattern;
}

/**
 * @brief Produce a random reassortement pattern.
 *
 * @returns reassortement pattern
 */
boost::dynamic_bitset<> haploid_clone::reassortment_pattern() {
	boost::dynamic_bitset<> rec_pattern;
	//the blocks of the bitset are to long for the rng, hence divide them by 4
	int bpblock = rec_pattern.bits_per_block/4;
	long unsigned int temp_rec_pattern;
	int bits_left=number_of_loci, still_to_append;
	//set the random bitset with bpblock at a time
	while(bits_left>=4*bpblock){
		temp_rec_pattern=gsl_rng_uniform_int(evo_generator, 1<<bpblock);
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<bpblock;
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<(2*bpblock);
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<(3*bpblock);
		//cout <<bpblock<<"  "<<temp_rec_pattern<<endl;
		rec_pattern.append(temp_rec_pattern);
		bits_left-=4*bpblock;
	}
	//set the remaining bits one by one.
	//TODO: THIS NEEDS SOME EXCEPTION HANDLING SINCE bits_left can be too large
	still_to_append=bits_left;
	if (bits_left>=3*bpblock){
		temp_rec_pattern=gsl_rng_uniform_int(evo_generator, 1<<bpblock);
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<bpblock;
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<(2*bpblock);
		bits_left-=3*bpblock;
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bits_left)<<(3*bpblock);
	}
	else if (bits_left>=2*bpblock){
		temp_rec_pattern=gsl_rng_uniform_int(evo_generator, 1<<bpblock);
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<bpblock;
		bits_left-=2*bpblock;
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bits_left)<<(2*bpblock);
	}
	else if (bits_left>=bpblock){
		temp_rec_pattern=gsl_rng_uniform_int(evo_generator, 1<<bpblock);
		bits_left-=bpblock;
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bits_left)<<bpblock;
	}else{
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bits_left);
	}

	while (still_to_append){
		still_to_append--;
		//cout <<temp_rec_pattern<<"  "<<(temp_rec_pattern&(1<<bits_left))<<endl;
		rec_pattern.push_back(((temp_rec_pattern&(1<<still_to_append))>0));
	}
	return rec_pattern;
}

/**
 * @brief Produce and store a random sample of the population for stochastic processes.
 *
 * @param size size of the sample
 *
 * Produce a random sample of genotypes (labeled by their clone of origin)
 * which is to be used for mutations and sampling. The vector random_sample
 * exists for efficiency reasons, for it is faster to call the random number
 * generator once than many times.
 *
 * Note: if you want to get random clones, please use random_clone().
 */
void haploid_clone::produce_random_sample(int size) {
	if (HP_VERBOSE) cerr<<"haploid_clone::produce_random_sample(int): size "<<size<<"...";

	random_sample.reserve((size+50)*1.1);
	random_sample.clear();
	int thechosen,o;
	double frac = 1.1*(size+50)/population_size;
	//loop over all clones and choose a poisson distributed number of genoytpes
	for (unsigned int c=0; c<current_pop->size(); c++){
		thechosen=gsl_ran_poisson(evo_generator, frac*(*current_pop)[c].clone_size);
		//make sure it is not larger than the clone itself.
		thechosen=((*current_pop)[c].clone_size<thechosen)?((*current_pop)[c].clone_size):thechosen;
		//add each of the chosen individually to the random_sample vector
		if (thechosen)	for (o=0; o<thechosen; o++) random_sample.push_back(c);
	}
	gsl_ran_shuffle(evo_generator, &random_sample[0], random_sample.size(), sizeof(int));
	//random_shuffle(random_sample.begin(), random_sample.end());
	if (HP_VERBOSE) cerr<<"done"<<endl;
}

/**
 * @brief Get a random clone from the population.
 *
 * @returns the index of the random clone.
 *
 * The probability density function from which the individual is chosen is flat over the
 * population (larger clones are proportionally more likely to be returned here).
 *
 * Note: for efficiency reasons, the class keeps a storage of random clone indices.
 * If you need much more than, say, 1000 random clones at a time, please call
 * produce_random_sample with the required size in advance.
 */
int haploid_clone::random_clone() {
	int rclone;
	int size = 1000;
	if (random_sample.size()>1){
		rclone=random_sample.back();
		random_sample.pop_back();
		return rclone;
	}else{	//if no genotypes left in sample, produce new sample
		produce_random_sample(min(population_size, size));
		rclone=random_sample.back();
		random_sample.pop_back();
		return rclone;
	}
}

/**
 * @brief Sample random indivduals from the population.
 *
 * @param n_o_individuals number of individuals to sample
 * @param sample pointer to vector where to put the result.
 *
 * @returns zero if successful, nonzero otherwise
 *
 * The results are not returned as a vector for performance reasons, as one might
 * want to get a lot of random clones. *sample may be not empty (but must be allocated).
 * In any case, clone numbers of the sampled individuals are appended to *sample. Hence,
 * you can use this function iteratively (although there might not be a good reason to
 * do so).
 */
int haploid_clone::random_clones(unsigned int n_o_individuals, vector <int> *sample) {
	sample->reserve(n_o_individuals);
	for(unsigned int i=0; i< n_o_individuals; i++)
		sample->push_back(random_clone());
	return 0;
}

/**
 * @brief Add the genotype specified by a bitset to the current population in in n copies
 *
 * @param genotype genotype being added
 * @param n number of copies of that genotype
 *
 * Note: this function also calculates the traits and fitness of the new individual.
 */
void haploid_clone::add_genotypes(boost::dynamic_bitset<> genotype, int n) {
	clone_t tempgt(number_of_traits);
	tempgt.genotype=genotype;
	tempgt.clone_size = n;
	calc_individual_fitness(&tempgt);
	current_pop->push_back(tempgt);
	population_size+=n;
}

/**
 * @brief Get the log of the exp-average fitness plus relaxation term.
 *
 * @returns the carrying capacity, i.e. the current mean fitness plus term that causes the population size to relax to carrying_capacity
 */
double haploid_clone::relaxation_value() {
	if (HP_VERBOSE) {cerr <<"haploid_clone::relaxation_value()..."<<endl;}

//	// FIXME: OLD
//	update_traits();
//	update_fitness();
//	calc_fitness_stat();
//	double relax = fitness_stat.mean + (MIN(0.6931*(double(population_size)/carrying_capacity-1),2.0));
//	if (HP_VERBOSE) {cerr<<"mean fitness = "<<fitness_stat.mean<<"..."<<endl;}

	double fitness_max = get_max_fitness();
	double logmean_expfitness = get_logmean_expfitness(fitness_max);
	// the second term is the growth rate when we start from N << carrying capacity
	double relax = logmean_expfitness + (fmin(0.6931*(double(population_size)/carrying_capacity-1),2.0)) + fitness_max;
	if (HP_VERBOSE) {cerr<<"log(<exp(F-Fmax)>) = "<<logmean_expfitness<<"..."<<endl;}

	if (HP_VERBOSE) {cerr<<"relaxation value = "<<relax<<"...done."<<endl;}
	return relax;
}

/**
 * @brief Get the fitness of the fittest individual.
 *
 * @returns the fitness of the fittest clone
 *
 * Note: this function recalculates the max fitness.
 */
double haploid_clone::get_max_fitness() {
	double fitness_max = (*current_pop)[0].fitness;
	for (unsigned int i=0; i<current_pop->size(); i++)
		fitness_max = fmax(fitness_max, (*current_pop)[i].fitness);
	return fitness_max;
}

/**
 * @brief Calculate and store fitness population statistics.
 *
 * The fitness statistics are stored in fitness_stat.
 *
 * Note: this function assumes that fitness is up to date. If you are not sure, call update_traits() and update_fitness() first.
 */
void haploid_clone::calc_fitness_stat() {
	if (HP_VERBOSE) {cerr <<"haploid_clone::calc_fitness_stat()...";}
	double temp,temp1;
	int t,t1, csize;
	fitness_stat.mean=0;
	fitness_stat.variance=0;
	population_size=0;
	//loop over clones and add stuff up
	for (unsigned int c=0; c<current_pop->size(); c++){
		csize = (*current_pop)[c].clone_size;
		temp=(*current_pop)[c].fitness;
		//if (HP_VERBOSE) {cerr <<"fitness["<<c<<"] = "<<temp<<"...";}
		fitness_stat.mean+=temp*csize;
		fitness_stat.variance+=temp*temp*csize;
		population_size+=csize;
	}
	if (population_size==0){
		cerr <<"haploid_clone::calc_fitness_stat(): population extinct! clones: "<<current_pop->size()<<endl;
	}
	//if (HP_VERBOSE) {cerr <<"pop size: "<<population_size<<", sum of fitnesses: "<<fitness_stat.mean<<"...";}
	fitness_stat.mean/=population_size;
	fitness_stat.variance/=population_size;
	fitness_stat.variance-=fitness_stat.mean*fitness_stat.mean;
	if (HP_VERBOSE) {cerr <<"done."<<endl;}
}


/**
 * @brief Get the population exp-average of fitness, which is used for keeping the population size fixed.
 *
 * @param fitness_max Maximal fitness in the population, used to avoid explosion of exponentials.
 * 
 * @returns the population exp-average of fitness
 *
 * Mathematically, this is \f$ \log\left( \left< e^(F-F_{max}) \right> \tight) \f$. The baseline is \f$ F_{max} \f$
 * in order to avoid exponentiating large numbers.
 *
 */
double haploid_clone::get_logmean_expfitness(double fitness_max) {
	if (HP_VERBOSE) {cerr <<"haploid_clone::get_logmean_expfitness()...";}
	double logmean_expfitness = 0;
	//loop over clones and add stuff up
	for (unsigned int c=0; c<current_pop->size(); c++){
		logmean_expfitness += (*current_pop)[c].clone_size * exp((*current_pop)[c].fitness - fitness_max);
	}
	logmean_expfitness /= population_size;
	logmean_expfitness = log(logmean_expfitness);
	if (HP_VERBOSE) {cerr <<"done."<<endl;}
	return logmean_expfitness;
}


/**
 * @brief Calculate and store trait population statistics and covariances
 *
 * The traits statistics are stored in traits_stat.
 *
 * Note: this function assumes that traits are up to date. If you are not sure, call update_traits() first.
 */
void haploid_clone::calc_trait_stat() {
	if (HP_VERBOSE) {cerr <<"haploid_clone::calc_trait_stat()...";}
	double temp,temp1;
	int t,t1, csize;
	population_size=0;

	// reset class attributes
	for(t=0; t<number_of_traits; t++){
		trait_stat[t].mean=0;
		trait_stat[t].variance=0;
		for(t1=0; t1<number_of_traits; t1++){
			trait_covariance[t][t1]=0;
		}
	}

	//loop over clones and add stuff up
	for (unsigned int c=0; c<current_pop->size(); c++){
		csize = (*current_pop)[c].clone_size;
		for(t=0; t<number_of_traits; t++){
			temp=(*current_pop)[c].trait[t];
			trait_stat[t].mean+=temp*csize;
			trait_stat[t].variance+=temp*temp*csize;
			for(t1=0; t1<number_of_traits; t1++){
				temp1=(*current_pop)[c].trait[t1];
				trait_covariance[t][t1]+=temp*temp1*csize;
			}
		}
		temp=(*current_pop)[c].fitness;
		population_size+=csize;
	}
	//complain if population went extinct
	if (population_size==0){
		cerr <<"haploid_clone::calc_trait_stat(): population extinct! clones: "
				<<current_pop->size()<<endl;
	}
	//normalize the means and variances
	for(t=0; t<number_of_traits; t++){
		trait_stat[t].mean/=population_size;
		trait_stat[t].variance/=population_size;
		trait_stat[t].variance -= trait_stat[t].mean*trait_stat[t].mean;
	}
	//calculate the covariances
	for(t=0; t<number_of_traits; t++){
		for(t1=0; t1<number_of_traits; t1++){
			trait_covariance[t][t1]/=population_size;
			trait_covariance[t][t1] -= trait_stat[t].mean*trait_stat[t1].mean;
		}
	}
	if (HP_VERBOSE) {cerr <<"done"<<endl;}
}

/**
 * @brief Print all allele frequency into a stream provided.
 *
 * @param out stream to put the allele frequencies (usually a file or stdout).
 *
 * @returns zero if successful, nonzero otherwise
 */
int haploid_clone::print_allele_frequencies(ostream &out) {
	if (out.bad())
	{
		cerr <<"haploid_clone::print_allele_frequencies: bad stream\n";
		return HP_BADARG;
	}
	calc_stat();
	out <<setw(10)<<generation;
	for (int l=0; l<number_of_loci; l++) out <<setw(15)<<allele_frequencies[l];
	out <<endl;
	return 0;
}

/**
 * @brief Read the output of Hudson's ms and use it to initialize the genotype distribution.
 *
 * @param gts genotypes output of _ms_
 * @param skip_locus positionof the locus to be skipped
 * @param multiplicity number of times each genotype must be added
 *
 * @returns zero if successful, error codes otherwise
 *
 * ms loci are fed into the genotype with a locus that is skipped. Each ms genotype is added multiple times.
 *
 */
int haploid_clone::read_ms_sample(istream &gts, int skip_locus, int multiplicity) {
	if (gts.bad()){
		cerr<<"haploid_clone::read_ms_sample(): bad stream!\n";
		return HP_BADARG;
	}
	//line buffer to read in the ms input
	char *line= new char [2*number_of_loci+5000];
	bool found_gt=false;
	string header;
	int count=0;
	int segsites, site, locus;
	segsites=0;

	//new genotype to be read in from ms
	boost::dynamic_bitset<> newgt(number_of_loci);
	//reset population
	current_pop->clear();
	random_sample.clear();
	population_size=0;
	if (mem){
		init_genotypes(0);
		population_size=0;
		//loop over each line of the ms out-put file
		while(gts.eof()==false){
			gts.get(line, 2*number_of_loci+5000);
			gts.get();
			//skip over empty lines
			while (gts.peek()=='\n'){gts.get();};

			//cout <<count<<"  "<<gt<<" "<<line<<endl;
			count++;
			//if the first genotype has been found
			if (found_gt and line[0]!='\0'){
				newgt.reset();
				//go over the line and assing loci, skip over the "skip_locus"
				for (site=0; site<segsites; site++){
					if (line[site]=='1'){
					  if (site<skip_locus){
					    locus=site;
					    newgt.set(locus);
					  }
					  else{
					    locus=site+1;
					    newgt.set(locus);
					  }
					}
				}
				add_genotypes(newgt, multiplicity);
			}else{ //reading the header
				header.assign(line);
				if (header.compare(0,2, "//")==0) {
					gts.get(line, 2*number_of_loci);
					gts.get();
					while (gts.peek()=='\n'){gts.get();};
					cerr <<count<<"  "<<found_gt<<" "<<line<<endl;
					header.assign(line);
					segsites=atoi(header.substr(9,header.size()-9).c_str());
					gts.get(line, 2*number_of_loci);
					gts.get();
					while (gts.peek()=='\n'){gts.get();};
					cerr <<count<<"  "<<found_gt<<" "<<line<<endl;
					found_gt=true;
				}
			}
		}
	}
	delete [] line;
	return 0;
}

/**
 *
 * @brief Read the output of Hudson's ms and use it to initialize the genotype distribution.
 *
 * @param gts genotypes output of _ms_
 * @param skip_locus positionof the locus to be skipped
 * @param multiplicity number of times each genotype must be added
 * @param distance distance between two polymorphic sites
 *
 * @returns zero if successful, error codes otherwise
 *
 * ms loci are fed into the genotype at distance "distance", i.e. there are distance-1 monomorphic loci.
 * One locus is skipped. Each ms genotype is added multiple times.
 */
int haploid_clone::read_ms_sample_sparse(istream &gts, int skip_locus, int multiplicity, int distance) {
	if (gts.bad()){
		cerr<<"haploid_clone::read_ms_sample(): bad stream!\n";
		return HP_BADARG;
	}
	//line buffer to read in the ms input
	char *line= new char [2*number_of_loci+5000];
	bool found_gt=false;
	string header;
	int count=0;
	int segsites, site, locus;
	segsites=0;

	//new genotype to be read in from ms
	boost::dynamic_bitset<> newgt(number_of_loci);
	//reset population
	current_pop->clear();
	random_sample.clear();
	population_size=0;
	if (mem){
		init_genotypes(0);
		population_size=0;
		//loop over each line of the ms out-put file
		while(gts.eof()==false){
			gts.get(line, 2*number_of_loci+5000);
			gts.get();
			//skip over empty lines
			while (gts.peek()=='\n'){gts.get();};

			//cout <<count<<"  "<<gt<<" "<<line<<endl;
			count++;
			//if the first genotype has been found
			if (found_gt and line[0]!='\0'){
				newgt.reset();
				//go over the line and assing loci, skip over the "skip_locus"
				for (site=0; site<segsites and site*distance<get_number_of_loci(); site++){
					if (line[site]=='1'){
					    locus=site*distance;
					    if (locus!=skip_locus){
						    newgt.set(locus);
					    }
					}
				}
				add_genotypes(newgt, multiplicity);
			}else{ //reading the header
				header.assign(line);
				if (header.compare(0,2, "//")==0) {
					gts.get(line, 2*number_of_loci);
					gts.get();
					while (gts.peek()=='\n'){gts.get();};
					cerr <<count<<"  "<<found_gt<<" "<<line<<endl;
					header.assign(line);
					segsites=atoi(header.substr(9,header.size()-9).c_str());
					gts.get(line, 2*number_of_loci);
					gts.get();
					while (gts.peek()=='\n'){gts.get();};
					cerr <<count<<"  "<<found_gt<<" "<<line<<endl;
					found_gt=true;
				}
			}
		}
	}
	delete [] line;
	return 0;
}


/**
 * @brief Calculate Hamming distance between two sequences.
 *
 * @param gt1 first sequence
 * @param gt2 second sequence
 * @param chunks (pointer to) vector of ranges (C pairs), e.g. ((0,10), (25, 28), (68, 70))
 * @param every check only every X sites, starting from the first of each chunk
 *
 * *Note*: you cannot use `every` without `chunks`.
 *
 * *Note*: every chunk is a C array of length 2, i.e. a pointer to lower boundary of the interval. Thus we have the following situation:
 * - `chunks` is of type vector `<unsigned int *> *`
 * - `*chunks` is of type vector `<unsigned int *>`
 * - `(*chunks)[i]` is of type `*unsigned int`
 * - `(*chunks)[i][0]` is of type `unsigned int`, and indicates the initial site for this chunk
 * - `(*chunks)[i][1]` is of type `unsigned int`, and indicates the (final site +1) for this chunk
 *
 * As a consequence, for each chunk, chunk[1] - chunk[0] is the length of the chunk.
 *
 * @returns Hamming distance, not normalized
 *
 * When you prepare the vector of chunks for this function, please use C++ methods (`new`) that take care of memory management, or be very careful about memory leaks.
 *
 * *Note*: this function is overloaded with simpler arguments (e.g. if you want to use clone indices).
 */
int haploid_clone::distance_Hamming(boost::dynamic_bitset<> gt1, boost::dynamic_bitset<> gt2, vector <unsigned int *> *chunks, unsigned int every) {
	// check whether we have chunks at all
	if((!chunks) || (chunks->size() == 0)) {
		if(every!=1) return HP_BADARG;
		else return (gt1 ^ gt2).count();
	}

	unsigned int d = 0;
	unsigned int pos;

	// check that the chunks make sense
	if((every<1) || (every>=number_of_loci)) return HP_BADARG;
	for(int i=0; i < chunks->size(); i++) {
		if((*chunks)[i][1] <= (*chunks)[i][0]) return HP_BADARG;
		if((*chunks)[i][1] >= number_of_loci) return HP_BADARG;
		if((*chunks)[i][0] >= number_of_loci) return HP_BADARG;
		
		for(pos = (*chunks)[i][0]; pos < (*chunks)[i][1]; pos+=every) {
			d += (unsigned int)(gt1[pos] != gt2[pos]);
		}
	}
	return d;
}


/**
 * @brief Calculate the cumulative partition of sequences into clones.
 *
 * @param partition_cum the vector to be filled
 *
 * @returns vector of cumulative clone sizes
 *
 * *Example*: if there are three clones of sizes (100, 22, 3) this function will
 * fill the vector (100, 122, 125). The last element is of course the population
 * size See also get_population_size.
 *
 * Note: the vector is taken in input by reference for performance reasons, since it can get huge.
 */
int haploid_clone::partition_cumulative(vector <unsigned int> &partition_cum) {	
	partition_cum.push_back((*current_pop)[0].clone_size);		
	for (size_t i = 1; i < get_number_of_clones(); i++) {
		partition_cum.push_back((*current_pop)[i].clone_size + partition_cum[i-1]);		
	}
	return 0;
}

/**
 * @brief Calculate mean and variance of the divergence from the [00...0] bitset.
 *
 * @param n_sample size of the statistical sample to use (the whole pop is often too large)
 *
 * @returns mean and variance of the divergence in a stat_t 
 */
stat_t haploid_clone::get_divergence_statistics(unsigned int n_sample)
{
	stat_t div;
	unsigned int tmp;
	vector <int> clones;
	produce_random_sample(n_sample);
	random_clones(n_sample, &clones);

	for (size_t i=0; i < n_sample; i++) {
		tmp = ((*current_pop)[clones[i]].genotype).count();
		div.mean += tmp;
		div.variance += tmp * tmp;
	}
	div.mean /= n_sample;
	div.variance /= n_sample;
	div.variance -= div.mean * div.mean;
	return div;
}

/**
 * @brief Calculate diversity in the current population (Hamming distance between pairs of sequences).
 *
 * @param n_sample size of the statistical sample to use (the whole pop is often too large)
 *
 * @returns mean and variance of the diversity in a stat_t
 */
stat_t haploid_clone::get_diversity_statistics(unsigned int n_sample)
{
	stat_t div;
	unsigned int tmp;
	vector <int> clones1;
	vector <int> clones2;
	produce_random_sample(n_sample * 2);
	random_clones(n_sample, &clones1);
	random_clones(n_sample, &clones2);

	for (size_t i=0; i < n_sample; i++) {
		if (clones1[i] != clones2[i]) {
			tmp = distance_Hamming(clones1[i],clones2[i]);
			div.mean += tmp;
			div.variance += tmp * tmp;
		}
	}
	div.mean /= n_sample;
	div.variance /= n_sample;
	div.variance -= div.mean * div.mean;
	return div;
}


/**
 * @brief Calculate histogram of fitness from traits.
 *
 * @param hist pointer to the gsl_histogram to fill
 * @param bins number of bins in the histogram 
 *
 * *Note*: the output histogram might have less bins than requested if the
 * sample size is too small.
 *
 * @param n_sample size of the random sample to use (the whole population is often too large)
 *
 *
 * @returns zero if successful, error codes otherwise
 *
 * There is a small problem here, namely that sometimes the fitness distribution has a
 * horrible tail which messes up the calculation of the bin width. Thus we first calculate
 * the fitness average and variance, and then set the bin width so that the deleterious
 * tail is within 2 sigma or so.
 *
 * *Note*: this function allocates memory for the histogram *only* if there is no binning error.
 * If you get HP_NOBINSERR, no memory is allocated. The user is expected to release the memory
 * manually.
 */
int haploid_clone::get_fitness_histogram(gsl_histogram **hist, unsigned int bins, unsigned int n_sample) {
	if (HP_VERBOSE) cerr <<"haploid_clone::get_fitness_histogram()...";

	// Calculate fitness of the sample
	double fitnesses[n_sample];
	vector <int> clones;
	produce_random_sample(n_sample);
	random_clones(n_sample, &clones);
	for(size_t i=0; i < n_sample; i++)
		fitnesses[i] = (*current_pop)[clones[i]].fitness;

	// Set the bins according to average and variance in fitness in the population
	calc_fitness_stat();
	double fitmean = fitness_stat.mean;
	double fitstd = sqrt(fitness_stat.variance);
	double histtail = 2 * fitstd;

	// Prepare histogram
	double fitmax = *max_element(fitnesses, fitnesses+n_sample);
	double fitmin = *min_element(fitnesses, fitnesses+n_sample);
	double fitmin_hist = fitmean - histtail;

	// Sometimes the population is homogeneous (neutral models)
	if(fitmin >= fitmax)
		return HP_NOBINSERR;

	//TODO: choose a decent criterion
	bins = min(n_sample / 30, bins);

	double width = (fitmax - fitmin_hist) / (bins - 1);
	*hist = gsl_histogram_alloc(bins); 
	gsl_histogram_set_ranges_uniform(*hist, fitmin_hist - 0.5 * width, fitmax + 0.5 * width);

	// Fill and scale histogram
	for (size_t i = 0; i < n_sample; i++)
		gsl_histogram_increment(*hist, fitnesses[i]);
	gsl_histogram_scale(*hist, 1/(double)n_sample);

	if (HP_VERBOSE) cerr <<"done"<<endl;

	// Of course, the user must take care of the memory allocated
	return 0;
}


/**
 * @brief Get histogram of divergence from the [00...0] bitset.
 *
 * @param hist pointer to the gsl_histogram to fill
 * @param bins number of bins in the histogram
 *
 * *Note*: an antialiasing algorithm is used to enforce integer bin widths.
 * *Note*: because of the antialiasing requirement, the last bin might be underrepresented.
 * *Note*: the number of used bins in the output might be smaller than requested, depending
 * on sample size and antialiasing.
 *
 * @param chunks (pointer to) vector of ranges (C pairs), e.g. ((0,10), (25, 28), (68, 70))
 * @param every check only every X sites, starting from the first of each chunk
 *
 * *Note*: you cannot use `every` without `chunks`.
 *
 * *Note*: every chunk is a C array of length 2, i.e. a pointer to lower boundary of the interval. Thus we have the following situation:
 * - `chunks` is of type vector `<unsigned int *> *`
 * - `*chunks` is of type vector `<unsigned int *>`
 * - `(*chunks)[i]` is of type `*unsigned int`
 * - `(*chunks)[i][0]` is of type `unsigned int`, and indicates the initial site for this chunk
 * - `(*chunks)[i][1]` is of type `unsigned int`, and indicates the (final site +1) for this chunk
 *
 * As a consequence, for each chunk, chunk[1] - chunk[0] is the length of the chunk.
 *
 * @param n_sample size of the random sample to use (the whole population is often too large)
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_clone::get_divergence_histogram(gsl_histogram **hist, unsigned int bins, vector <unsigned int *> *chunks, unsigned int every, unsigned int n_sample)
{
	if (HP_VERBOSE) {cerr <<"haploid_clone::get_divergence_histogram(gsl_histogram **hist, unsigned int bins, vector <unsigned int *> *chunks, unsigned int every, unsigned int n_sample)...";}

	boost::dynamic_bitset<> gt_wt(number_of_loci);	// the [00...0] bitset
	int temp;
	unsigned int divs[n_sample];
	vector <int> clones;
	produce_random_sample(n_sample);
	random_clones(n_sample, &clones);
	for(size_t i=0; i < n_sample; i++) {
		temp = distance_Hamming(gt_wt, (*current_pop)[clones[i]].genotype, chunks, every);
		// negative distances are error codes
		if(temp < 0) return temp;
		else divs[i] = temp;
	}

	// Prepare the histogram
	unsigned long dmax = *max_element(divs, divs + n_sample);
	unsigned long dmin = *min_element(divs, divs + n_sample);

	// Antialiasing
	int width, binsnew;
	if (dmin == dmax)
		width = 1;
	else {
		width = (dmax - dmin) / (bins-1);
		width += ((dmax - dmin)%(bins-1))?1:0;
	}
	binsnew = ((dmax - dmin) / width) + 1;
	if (binsnew > bins) {
		if (HP_VERBOSE) cerr<<"wrong bins!: "<<"binsnew: "<<binsnew<<", bins: "<<bins<<", delta: "<<(dmax - dmin)<<endl;
		return HP_WRONGBINSERR;
	}

	// Fill and scale histogram
	*hist = gsl_histogram_alloc(binsnew); 
	gsl_histogram_set_ranges_uniform(*hist, dmin - 0.5 * width, dmax + 0.5 * width);
	for (size_t i = 0; i < n_sample; i++)
		gsl_histogram_increment(*hist, divs[i]);
	gsl_histogram_scale(*hist, 1/(double)n_sample);
	
	if (HP_VERBOSE) cerr<<"done.";
	return 0;
}

/**
 * @brief Get histogram of diversity in the population (mutual Hamming distance).
 *
 * @param hist pointer to the gsl_histogram to fill
 * @param bins number of bins in the histogram
 *
 * *Note*: an antialiasing algorithm is used to enforce integer bin widths.
 * *Note*: because of the antialiasing requirement, the last bin might be underrepresented.
 * *Note*: the number of used bins in the output might be smaller than requested, depending
 * on sample size and antialiasing.
 *
 * @param chunks (pointer to) vector of ranges (C pairs), e.g. ((0,10), (25, 28), (68, 70))
 * @param every check only every X sites, starting from the first of each chunk
 *
 * *Note*: you cannot use `every` without `chunks`.
 *
 * *Note*: every chunk is a C array of length 2, i.e. a pointer to lower boundary of the interval. Thus we have the following situation:
 * - `chunks` is of type vector `<unsigned int *> *`
 * - `*chunks` is of type vector `<unsigned int *>`
 * - `(*chunks)[i]` is of type `*unsigned int`
 * - `(*chunks)[i][0]` is of type `unsigned int`, and indicates the initial site for this chunk
 * - `(*chunks)[i][1]` is of type `unsigned int`, and indicates the (final site +1) for this chunk
 *
 * As a consequence, for each chunk, chunk[1] - chunk[0] is the length of the chunk.
 *
 * @param n_sample size of the random sample to use (the whole population is often too large)
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_clone::get_diversity_histogram(gsl_histogram **hist, unsigned int bins, vector <unsigned int *> *chunks, unsigned int every, unsigned int n_sample)
{
	if (HP_VERBOSE) {cerr <<"haploid_clone::get_diversity_histogram(gsl_histogram **hist, unsigned int bins, vector <unsigned int *> *chunks, unsigned int every, unsigned int n_sample)...";}
	int temp;
	unsigned int divs[n_sample];
	vector <int> clones1;
	vector <int> clones2;
	produce_random_sample(n_sample * 2);
	random_clones(n_sample, &clones1);
	random_clones(n_sample, &clones2);
	for(size_t i=0; i < n_sample; i++) {
		temp = distance_Hamming(clones1[i], clones2[i], chunks, every);
		// negative distances are error codes
		if(temp < 0) return temp;
		else divs[i] = temp;
	}

	// Prepare the histogram
	unsigned long dmax = *max_element(divs, divs + n_sample);
	unsigned long dmin = *min_element(divs, divs + n_sample);

	// Aliasing
	int width, binsnew;
	if (dmin == dmax)
		width = 1;
	else {
		width = (dmax - dmin) / (bins -1);
		width += ((dmax - dmin)%(bins-1))?1:0;
	}
	binsnew = ((dmax - dmin) / width) + 1;
	if (binsnew > bins) {
		if (HP_VERBOSE) cerr<<"wrong bins!: "<<"binsnew: "<<binsnew<<", bins: "<<bins<<", delta: "<<(dmax - dmin)<<endl;
		return HP_WRONGBINSERR;
	}
	
	// Fill and scale histogram
	*hist = gsl_histogram_alloc(binsnew); 
	gsl_histogram_set_ranges_uniform(*hist, dmin - 0.5 * width, dmax + 0.5 * width);
	for (size_t i = 0; i < n_sample; i++)
		gsl_histogram_increment(*hist, divs[i]);
	gsl_histogram_scale(*hist, 1/(double)n_sample);
	
	if (HP_VERBOSE) cerr<<"done.";
	return 0;
}
