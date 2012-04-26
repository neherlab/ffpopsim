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
 * The sequence is assumed to be linear (not circular).
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
 * Memory is released here. This needs to be virtual, because subclasses might want to free their own memory
 * (but is a small detail in our typical usage cases). Subclasses invoke the baseclass destructor at the end anyway
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
	target_pop_size=N_in;

	//In case no seed is provided use current second and add process ID
	if (rng_seed==0){
		seed=time(NULL)+getpid();
	}else{seed=rng_seed;}

	mem=false;
	int err=allocate_mem();
	return err;
}

/*
 * Allocate all the necessary memory, initialze the RNG
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

// create genotypes at random from allele frequencies saved in nu[k] k=0..number_of_loci-1
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
	// if number of genotypes to be drawn is not specified, use default target population size
	if (n_o_genotypes==0) ngt=target_pop_size;
	else ngt=n_o_genotypes;

	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes(double* nu, int n_o_genotypes) reset current\n";
	current_pop->clear();	//reset the current population
	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes(double* nu, int n_o_genotypes) reset random sample\n";
	random_sample.clear();	//and the random sample
	boost::dynamic_bitset<> tempgt(number_of_loci);
	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes(double* nu, int n_o_genotypes) add "<<ngt<<" genotypes of length "<<number_of_loci<<"\n";
	for (i=0; i<ngt; i++)
	{
		tempgt.reset();
		for(locus=0; locus<number_of_loci; locus++){	//set all loci
			if (gsl_rng_uniform(evo_generator)<nu[locus])	tempgt.set(locus);
		}
		add_genotypes(tempgt,1);	//add genotype with multiplicity 1
		if (HP_VERBOSE) cerr <<i<<endl;
	}

	//calculate its fitness and recombination rates
	calc_stat();
	if (HP_VERBOSE) cerr <<"done.\n";
	//set generations counter to zero and calculate the statistics for the intial configuration
	generation=0;
	return 0;
}



/**
 * @brief Initialize with a single genotypic clone (00...0)
 *
 * @param n_o_genotypes size of the clone. If not chosen, use target_pop_size.
 *
 * @returns 0 if successful, nonzero otherwise.
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
	pop_size=0;
	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes(int n_o_genotypes) with population size "<<n_o_genotypes<<"...";
	generation=0;
	current_pop->clear();
	random_sample.clear();
	boost::dynamic_bitset<> tempgt(number_of_loci);
	tempgt.reset();					//set all bits to zero
	if (n_o_genotypes>=0){				//add n_o_genotypes copies of this genotype
		add_genotypes(tempgt,n_o_genotypes);
	}else{
		add_genotypes(tempgt,target_pop_size);
	}

	calc_stat();					//make sure everything is calculated
	if (HP_VERBOSE) cerr <<"done."<<endl;
	return 0;
}


/**
 * @brief calculate allele frequencies and mean and variance in fitness
 */
void haploid_clone::calc_allele_freqs()
{
	if (HP_VERBOSE) cerr<<"haploid_clone::calc_allele_freqs()...";
	int locus,cs;
	pop_size=0;
	for (locus=0; locus<number_of_loci; allele_frequencies[locus++] = 0);
	//loop over all clones
	for (unsigned int i=0; i<current_pop->size(); i++)
	{
		cs=(*current_pop)[i].clone_size;
		for (locus=0; locus<number_of_loci; locus++)	//add clone size to allele frequency of clone carries allele
			if ((*current_pop)[i].genotype[locus])
				allele_frequencies[locus] += cs;
		pop_size += cs;
	}
	//convert counts into frequencies
	for (locus=0; locus<number_of_loci; allele_frequencies[locus++] /= pop_size);
	if (HP_VERBOSE) cerr<<"done.\n";
}

double haploid_clone::get_pair_frequency(int locus1, int locus2)
{
	if (HP_VERBOSE) cerr<<"haploid_clone::get_pair_frequency()...";
	double frequency = 0;
	for (unsigned i=0; i<current_pop->size(); i++)
		if ((*current_pop)[i].genotype[locus1] and (*current_pop)[i].genotype[locus2])
		        frequency += (*current_pop)[i].clone_size;
	frequency /= pop_size;
	if (HP_VERBOSE) cerr<<"done.\n";
	return frequency;
}

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
 * @brief evolve for g generations
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
 */
int haploid_clone::evolve(int gen){
	int err=0, gtemp=0;
	if (HP_VERBOSE) cerr<<"haploid_clone::evolve(int gen)...";
	while((err == 0) && (gtemp < gen)) {
		if (HP_VERBOSE) cerr<<"generation "<<(generation+gtemp)<<endl;
		random_sample.clear();			//discard the old random sample
		if(err==0) err=select_gametes();	//select a new set of gametes (partitioned into sex and asex)
		if(err==0) err=add_recombinants();	//do the recombination between pairs of sex gametes
		if(err==0) err=swap_populations();	//make the new population the current population
		if(err==0) err=mutate();		//mutation step
		gtemp++;
	}
	generation+=gtemp;
	if (HP_VERBOSE) {
		if(err==0) cerr<<"done."<<endl;
		else cerr<<"error "<<err<<"."<<endl;
	}
	return err;
}

/**
 * @brief selection step
 *
 * Random Poisson offspring numbers are drawn from all parents, proportionally to their fitness.
 */
int haploid_clone::select_gametes()
{
	if (HP_VERBOSE) cerr<<"haploid_clone::select_gametes()...";
	//determine the current mean fitness, which includes a term to keep the population size constant
	double cpot=chemical_potential();
	//draw gametes according to parental fitness
	int os,o, nrec;
	int err=0;
	pop_size=0;
	
	//to speed things up, reserve the expected amount of memory for sex gametes and the new population (+10%)
	sex_gametes.reserve(pop_size*outcrossing_probability*1.1);
	sex_gametes.clear();
	new_pop->reserve(current_pop->size()*1.1);
	new_pop->clear();
	for (unsigned int i=0; i<current_pop->size(); i++) {
		//if (HP_VERBOSE) cerr<<i<<" ";
		//poisson distributed random numbers -- mean exp(f)/bar{exp(f)}) (since death rate is one, the growth rate is (f-bar{f})
		if ((*current_pop)[i].clone_size>0){
			//the number of asex offspring of clone[i] is poisson distributed around e^{F-mF}(1-r)
			//if (HP_VERBOSE) cerr<<i<<": relative fitness = "<<((*current_pop)[i].fitness-cpot)<<", Poisson intensity = "<<((*current_pop)[i].clone_size*exp((*current_pop)[i].fitness-cpot)*(1-outcrossing_probability))<<endl;
			os=gsl_ran_poisson(evo_generator, (*current_pop)[i].clone_size*exp((*current_pop)[i].fitness-cpot)*(1-outcrossing_probability));

			if (os>0){
				// clone[i] to new_pop with os as clone size
				new_pop->push_back((*current_pop)[i]);
				new_pop->back().clone_size=os;
				pop_size+=os;
			}
			//draw the number of sexual offspring, add them to the list of sex_gametes one by one
			if (outcrossing_probability>0){
				nrec=gsl_ran_poisson(evo_generator, (*current_pop)[i].clone_size*exp((*current_pop)[i].fitness-cpot)*outcrossing_probability);
				for (o=0; o<nrec; o++) sex_gametes.push_back(i);
			}
		}
	}
	if(pop_size<1) err = HP_EXTINCTERR;
	if (HP_VERBOSE) {
		if(err==0) cerr<<"done."<<endl;
		else {
			cerr<<"error "<<err<<".";
			if(err == HP_EXTINCTERR) cerr<<" The population went extinct!";
			cerr<<endl;
		}
	}
	return err;
}

/**
 * @brief mutation step (all loci)
 *
 * For efficiency reasons, if the mutation rate is zero, an if cycle skips the body altogether.
 * The user can therefore call this function safely even if in non-mutating populations, without
 * loss of performance.
 */
int haploid_clone::mutate() {
	if (HP_VERBOSE)	cerr <<"haploid_clone::mutate() ..."<<endl;
	int i, actual_n_o_mutations,locus=0;
	if(mutation_rate) {
		produce_random_sample(number_of_loci*mutation_rate*pop_size*2);
		for (locus = 0; locus<number_of_loci; locus++) {
			//draw the number of mutation that are to happen at this locus
			actual_n_o_mutations = gsl_ran_poisson(evo_generator, mutation_rate*pop_size);
			//introduce these mutations one by one
			for (i = 0; i < actual_n_o_mutations; ++i) {
				flip_single_locus(random_clone(), locus);
			}
		}
	}
	if (HP_VERBOSE)	cerr <<"done."<<endl;;
	return 0;
}

/**
 * @brief flip a spin at locus in random individual
 *
 * @param locus locus to flip
 *
 * @returns the (random) flipped individual clone
 */
int haploid_clone::flip_single_locus(int locus)
{
	int clonenum = random_clone();
	flip_single_locus(clonenum, locus);
	return clonenum;
}


//flip a spin (locus) in individual, assign new fitness and recombination rate.
void haploid_clone::flip_single_locus(unsigned int clonenum, int locus)
{
	//produce new genotype
	clone_t tempgt(number_of_traits);
	tempgt.genotype.resize(number_of_loci);
	tempgt.genotype= (*current_pop)[clonenum].genotype;
	//new clone size == 1, old clone reduced by 1 TODO:check that this makes sense
	tempgt.clone_size=1;
	(*current_pop)[clonenum].clone_size--;
	//flip the locus in new clone and calculate fitness
	tempgt.genotype.flip(locus);
	calc_individual_fitness(&tempgt);
	//add clone to current population
	current_pop->push_back(tempgt);
//FIXME	if (HP_VERBOSE) cerr <<"subpop::flip_single_spin(): mutated individual in clone "<<clonenum<<" at locus "<<locus<<endl;
}

/**
 * @brief pair and mate sexual gametes.
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
			//in new_genotypes[new_pop_size] and [new_pop_size+1]
			n_o_c=recombine(parent1, parent2);
		}
	}
	return 0;
}

/**
 * @brief make the the temporary population the current one
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
int haploid_clone::recombine(int parent1, int parent2)
{
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
	//cout <<rec_pattern<<endl;
	//cout <<(*current_pop)[parent1].genotype<<endl;
	//cout <<(*current_pop)[parent2].genotype<<endl;
	//cout <<offspring1.genotype<<endl;
	//cout <<offspring2.genotype<<endl<<endl;

	calc_individual_fitness(&offspring1);
	calc_individual_fitness(&offspring2);

	//add genotypes to new population
	new_pop->push_back(offspring1);
	new_pop->push_back(offspring2);
	pop_size+=2;
	return 1;
}

/**
 * @brief for each clone, update fitness assuming traits are already up to date
 */
void haploid_clone::update_fitness() {
	for (unsigned int i=0; i<current_pop->size(); i++) {
		calc_individual_fitness_from_traits(&(*current_pop)[i]);
	}
}

/**
 * @brief for each clone, recalculate its traits
 */
void haploid_clone::update_traits(){
	for (unsigned int i=0; i<current_pop->size(); i++) {
		calc_individual_traits(&(*current_pop)[i]);
	}
}

/**
 * @brief calculate trait and fitness statistics and allele frequences
 *
 * Four things are done in a row:
 * 1. traits are updated
 * 2. fitness is updated based on the new traits
 * 3. statistics of traits and fitness are calculated
 * 4. allele freqs are calculated
 * This function is therefore quite expensive. Please use its subblocks
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
 * @brief Calculate traits of the chosen clone
 *
 * @param tempgt clone whose traits are to be calculated
 */
void haploid_clone::calc_individual_traits(clone_t *tempgt){
	for (int t=0; t<number_of_traits; t++){
		tempgt->trait[t] = trait[t].get_func(&(tempgt->genotype));
	}
}

/**
 * @brief calculate fitness of a particular clone
 *
 * @param tempgt clone whose fitness is being calculated
 *
 * Note: this function also updates the traits information for the same clone,
 * because the phenotype is needed to calculate fitness.
 */
void haploid_clone::calc_individual_fitness(clone_t *tempgt){
	//calculate the new fitness value of the mutant
	calc_individual_traits(tempgt);
	calc_individual_fitness_from_traits(tempgt);
}

/*
 * choose a certain number of crossover points and produce a 0000111100011111101101 crossover pattern
 */
boost::dynamic_bitset<> haploid_clone::crossover_pattern(){
	int n_o_c=0;

	double total_rec = number_of_loci*crossover_rate;
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

/*
 * produce a random reassortement pattern
 */
boost::dynamic_bitset<> haploid_clone::reassortment_pattern(){
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
	//TODO
	//THIS NEEDS SOME EXCEPTION HANDLING SINCE bits_left can be too large
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
 * @brief Sample the population for stochastic processes.
 *
 * @param size size of the sample
 *
 * produce a random sample of genotypes (labeled by their clone of origin)
 * which is to be used for mutations and sampling. The vector random_sample
 * exists for efficiency reasons, for it is faster to call the random number
 * generator once than many times.
 *
 */
void haploid_clone::produce_random_sample(int size){
	if (HP_VERBOSE) cerr<<"haploid_clone::produce_random_sample(int): size "<<size<<"...";

	random_sample.reserve((size+50)*1.1);
	random_sample.clear();
	int thechosen,o;
	double frac = 1.1*(size+50)/pop_size;
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
 * @brief get a random clone from the population
 *
 * @returns 
 *
 * return the clone index of a random individual.
 *
 * The probability density function from which the individual is chosen is flat over the
 * population (larger clones are proportionally more likely to be returned here).
 *
 * The index is then removed from the random_sample vector, which, for efficiency reasons,
 * keeps a list of random indices ready for use.
 */
int haploid_clone::random_clone() {
	int rclone;
	int size = 1000;
	if (random_sample.size()>1){
		rclone=random_sample.back();
		random_sample.pop_back();
		return rclone;
	}else{	//if no genotypes left in sample, produce new sample
		produce_random_sample(min(pop_size, size));
		rclone=random_sample.back();
		random_sample.pop_back();
		return rclone;
	}
}

/**
 * @brief Sample random indivduals from the population
 *
 * @param n_o_individuals number of individuals to sample
 * @param sample pointer to vector where to put the result.
 *
 * @returns zero if successful, nonzero otherwise
 *
 * The results are not returned as a vector to avoid copying (for performance reasons),
 * and *sample may be not empty (but must be allocated). In any case, clone numbers of
 * the sampled individuals are appended to *sample. Hence, you can use this function
 * iteratively (although there might not be a good reason to do so).
 */
int haploid_clone::random_clones(unsigned int n_o_individuals, vector <int> *sample) {
	sample->reserve(n_o_individuals);
	for(unsigned int i=0; i< n_o_individuals; i++)
		sample->push_back(random_clone());
	return 0;
}

/*
 * add the genotype specified by a bitset to the current population in in n copies
 */
void haploid_clone::add_genotypes(boost::dynamic_bitset<> genotype, int n) {
	clone_t tempgt(number_of_traits);
	tempgt.genotype=genotype;
	tempgt.clone_size = n;
	calc_individual_fitness(&tempgt);
	current_pop->push_back(tempgt);
	pop_size+=n;
}

/**
 * @brief get the mean fitness plus relaxation term.
 *
 * @return the chemical potential, i.e. the current mean fitness plus term that causes the population size to relax to target_pop_size
 */
double haploid_clone::chemical_potential() {
	if (HP_VERBOSE) {cerr <<"haploid_clone::chemical_potential()..."<<endl;}
	double chem_pot;
	// Note: the update_X stuff is required, because the chemical potential is a global quantity
	update_traits();
	update_fitness();
	calc_fitness_stat();
	chem_pot = fitness_stat.mean + (MIN(0.6931*(double(pop_size)/target_pop_size-1),2.0));
	if (HP_VERBOSE) {cerr<<"mean fitness = "<<fitness_stat.mean<<"..."<<endl;}
	if (HP_VERBOSE) {cerr<<"cpot = "<<chem_pot<<"...done."<<endl;}
	return chem_pot;
}


/**
 * @brief get the fitness of the fittest individual
 *
 * @returns the fitness of the fittest clone
 */
double haploid_clone::get_max_fitness() {
	double mf=(*current_pop)[0].fitness;
	for (unsigned int i=0; i<current_pop->size(); i++)
		mf = MAX(mf, (*current_pop)[i].fitness);
	return mf;
}


/**
 * @brief calculate fitness population statistics
 *
 * Note: this function assumes that fitness is up to date. If you are not sure, call update_traits() and update_fitness() first.
 */
void haploid_clone::calc_fitness_stat() {
	if (HP_VERBOSE) {cerr <<"haploid_clone::calc_fitness_stat()...";}
	double temp,temp1;
	int t,t1, csize;
	fitness_stat.mean=0;
	fitness_stat.variance=0;
	pop_size=0;
	//loop over clones and add stuff up
	for (unsigned int c=0; c<current_pop->size(); c++){
		csize = (*current_pop)[c].clone_size;
		temp=(*current_pop)[c].fitness;
		//if (HP_VERBOSE) {cerr <<"fitness["<<c<<"] = "<<temp<<"...";}
		fitness_stat.mean+=temp*csize;
		fitness_stat.variance+=temp*temp*csize;
		pop_size+=csize;
	}
	if (pop_size==0){
		cerr <<"haploid_clone::calc_fitness_stat(): population extinct! clones: "<<current_pop->size()<<endl;
	}
	//if (HP_VERBOSE) {cerr <<"pop size: "<<pop_size<<", sum of fitnesses: "<<fitness_stat.mean<<"...";}
	fitness_stat.mean/=pop_size;
	fitness_stat.variance/=pop_size;
	fitness_stat.variance-=fitness_stat.mean*fitness_stat.mean;
	if (HP_VERBOSE) {cerr <<"done."<<endl;}
}

/**
 * @brief calculate trait population statistics and covariances
 *
 * Note: this function assumes that traits are up to date. If you are not sure, call update_traits() first.
 */
void haploid_clone::calc_trait_stat() {
	if (HP_VERBOSE) {cerr <<"haploid_clone::calc_trait_stat()...";}
	double temp,temp1;
	int t,t1, csize;
	pop_size=0;

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
		pop_size+=csize;
	}
	//complain if population went extinct
	if (pop_size==0){
		cerr <<"haploid_clone::calc_trait_stat(): population extinct! clones: "
				<<current_pop->size()<<endl;
	}
	//normalize the means and variances
	for(t=0; t<number_of_traits; t++){
		trait_stat[t].mean/=pop_size;
		trait_stat[t].variance/=pop_size;
		trait_stat[t].variance -= trait_stat[t].mean*trait_stat[t].mean;
	}
	//calculate the covariances
	for(t=0; t<number_of_traits; t++){
		for(t1=0; t1<number_of_traits; t1++){
			trait_covariance[t][t1]/=pop_size;
			trait_covariance[t][t1] -= trait_stat[t].mean*trait_stat[t1].mean;
		}
	}
	if (HP_VERBOSE) {cerr <<"done"<<endl;}
}

/*
 * convenience function, printing all allele frequency into a stream provided
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

/*
 * function that reads the output of Hudson's ms and uses it to initialize the genotype distribution.
 * ms loci are fed into the genotype with a locus that is skipped. Each ms genotype is added multiplicity times
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
	pop_size=0;
	if (mem){
		init_genotypes(0);
		pop_size=0;
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

/*
 * function that reads the output of Hudson's ms and uses it to initialize the genotype distribution.
 * ms loci are fed into the genotype at distance "distance", i.e. there are distance-1 monomorphic loci
 * one locus is skipped. Each ms genotype is added multiplicity times
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
	pop_size=0;
	if (mem){
		init_genotypes(0);
		pop_size=0;
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
 * @brief calculate Hamming distance between two sequences
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
 * @brief calculate the cumulative partition of sequences in the clones
 *
 * @returns vector of cumulative clone sizes
 *
 * *Example*: if there are three clones of sizes (100, 22, 3) this function will
 * return the vector (100, 122, 125). The last element is of course the population
 * size See also haploid_clone::get_pop_size().
 */
vector <unsigned int> haploid_clone::partition_cumulative()
{
	vector <unsigned int> partition_cum;
	partition_cum.push_back((*current_pop)[0].clone_size);		
	for (size_t i = 1; i < get_number_of_clones(); i++) {
		partition_cum.push_back((*current_pop)[i].clone_size + partition_cum[i-1]);		
	}
	return partition_cum;
}

/**
 * @brief Calculate mean and variance of the divergence from the [00...0] bitset
 *
 * @param n_sample size of the statistical sample to use (the whole pop is often too large)
 *
 * @returns mean and variance of the divergence in a struct_t 
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
 * @brief calculate diversity in the current population (Hamming distance between pairs of sequences)
 *
 * @param n_sample size of the statistical sample to use (the whole pop is often too large)
 *
 * @returns mean and variance of the diversity in a struct_t
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
 * @brief calculate histogram of fitness from traits
 *
 * @param hist pointer to the histogram to fill
 * @param bins number of bins in the histogram 
 * @param n_sample size of the random sample to use (the whole population is
 * often too large)
 *
 * *Note*: the output histogram might have less bins than requested if the
 * sample size is too small.
 *
 * @returns zero if successful, nonzero otherwise
 *
 * There is a small problem here, namely that the fitness distribution has a horrible tail which
 * messes up the calculation of the bin width. Thus we first calculate the fitness average and
 * variance, and then set the bin width so that the deleterious tail is within 2 sigma or so.
 *
 * *Note*: this function allocates memory for the histogram. The user is expected to release it
 * manually.
 */
int haploid_clone::get_fitness_histogram(gsl_histogram **hist, unsigned int bins, unsigned int n_sample) {
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
	double fmax = *max_element(fitnesses, fitnesses+n_sample);
	double fmin = fitmean - histtail;

	// Sometimes the population is homogeneous (neutral models)
	if(fmin >= fmax)
		return HP_NOBINSERR;

	bins = MIN(n_sample / 30, bins);	//TODO: choose a decent criterion
	double width = (fmax - fmin) / (bins - 1);
	*hist = gsl_histogram_alloc(bins); 
	gsl_histogram_set_ranges_uniform(*hist, fmin - 0.5 * width, fmax + 0.5 * width);

	// Fill and scale histogram
	for (size_t i = 0; i < n_sample; i++)
		gsl_histogram_increment(*hist, fitnesses[i]);
	gsl_histogram_scale(*hist, 1/(double)n_sample);

	// Of course, the user must take care of the memory allocated
	return 0;
}


/*
 * Get histogram of divergence, to ingestivate more in detail than just mean and variance
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

	// Aliasing
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

/*
 * Get histogram of diversity, to ingestivate more in detail than just mean and variance
 *
 * *Note*: the calculation of bin width and number requires an anti-aliasing algorithm.
 * This induces a (negative) bias in the last bin, but we cannot do anything for it (?).
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



