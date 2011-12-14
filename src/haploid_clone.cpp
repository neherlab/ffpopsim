// vim: tabstop=8:softtabstop=8:shiftwidth=8:noexpandtab
/*
 * haploid_clone.cpp
 *
 *  Created on: Aug 28, 2008
 *      Author: neher
 *    Modified: Dec 13, 2011
 *   Committer: Fabio Zanini
 */

#include "popgen.h"

haploid_clone::haploid_clone() {
	current_pop = &pop1;
	new_pop = &pop2;
	mem=false;
	cumulants_mem=false;
	circular=false;
}

haploid_clone::~haploid_clone() {
	if (mem) free_mem();
}

/**
 * Set up the population by specifying the population size, the length of the genome,
 * the random number generator seed and the number of traits for which genotype-trait maps
 * are allocated
 */
int haploid_clone::set_up(int N_in,int L_in,  int rng_seed, int n_o_traits)
{
	boost::dynamic_bitset<> temp;
	if (N_in<1 or L_in <1)
	{
		cerr <<"haploid_clone::set_up(): Bad Arguments!\n";
		return HP_BADARG;
	}

	if (8*sizeof(long int)!=temp.bits_per_block) {
		cerr <<"haploid_clone requires sizeof(long int) to be equal to bits_per_block of boost::dynamic_bitset";
		return HP_BADARG;
	}

	number_of_traits=n_o_traits;
	number_of_loci=L_in;
	number_of_individuals=2*N_in;
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

	if (HP_VERBOSE) cerr <<"haploid_clone::allocate_mem(): genotypes: "<<number_of_individuals<<"  loci: "<<number_of_loci<<endl;

	//Random number generator
	evo_generator=gsl_rng_alloc(RNG);
	gsl_rng_set(evo_generator, seed);
	//allocate all the memory
	genome = new int [number_of_loci+1];		// aux array holding range(0,number_of_loci) used to draw crossover points
	for (int i=0; i<number_of_loci; i++) genome[i]=i;
	crossovers= new int [number_of_loci];		// aux array holding crossover points

	if (HP_VERBOSE) cerr <<"allele frequencies...";
	allele_frequencies =new double [number_of_loci];
	gamete_allele_frequencies =new double [number_of_loci]; //allele frequencies after selection

	trait = new hypercube_function [number_of_traits];		//genotype trait function
	trait_stat = new stat [number_of_traits];				//structure holding trait statistics
	trait_covariance = new double* [number_of_traits];
	//initialize trait functions
	for (int t=0; t<number_of_traits;t++){
		trait[t].set_up(number_of_loci, gsl_rng_uniform_int(evo_generator, 1<<20));
		trait_covariance[t]=new double [number_of_traits];
	}

	mem=true;	//set memory flag to true
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
	if (n_o_genotypes<0 or n_o_genotypes>=number_of_individuals)
	{
		cerr <<"haploid_clone::init_genotypes(): number of genotypes has to be positive and smaller than max! got: "<<n_o_genotypes<<"!\n";
		return HP_BADARG;
	}
	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes(double* nu, int n_o_genotypes) ....";
	generation=0;

	// if number of genotypes to be drawn is not specified, use default target population size
	if (n_o_genotypes==0) pop_size=target_pop_size;
	else pop_size=n_o_genotypes;

	current_pop->clear();	//reset the current population
	random_sample.clear();	//and the random sample
	boost::dynamic_bitset<> tempgt(number_of_loci);
	for (i=0; i<pop_size; i++)
	{
		tempgt.reset();
		for(locus=0; locus<number_of_loci; locus++){	//set all loci
			if (gsl_rng_uniform(evo_generator)<nu[locus])	tempgt.set(locus);
		}
		add_genotypes(tempgt,1);	//add genotype with multiplicity 1
	}

	//calculate its fitness and recombination rates
	calc_fit();
	calc_stat();
	if (HP_VERBOSE) cerr <<"done.\n";
	//set generations counter to zero and calculate the statistics for the intial configuration
	generation=0;
	return 0;
}


//init with all allele frequencies equal to 0.5 -- n copies of each inital draw-default is 1
/*
int haploid_clone::init_genotypes_diverse(int n_o_genotypes, int no_copies)
{
	if (!mem)
	{
		cerr <<"haploid_clone::init_genotypes_diverse(): Allocate and set up first!\n";
		return HP_MEMERR;
	}
	if (n_o_genotypes<0 or n_o_genotypes>number_of_individuals)
	{
		cerr <<"haploid_clone::init_genotypes_diverse(): number of genotypes has to be positive and smaller than max! got: "<<n_o_genotypes<<"!\n";
		return HP_BADARG;
	}
	pop_size=0;
	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes_diverse(int n_o_genotypes, int no_copies) with population size "<<n_o_genotypes<<"...";
	generation=0;
	int i,j,w;
	for (i=0; i<n_o_genotypes;)
	{
		for(j=0;j<no_copies; j++);
		{
			for (w=0; w<number_of_words; w++)
			{		//random integer in 0:(1<<number_of_bits[k])-1 for each word
				genotypes[index(i,w)]=gsl_rng_uniform_int(evo_generator,(1<<number_of_bits[w]));
			}
			i++; pop_size++;
		}
	}
	fitness_values->number_of_values=pop_size;
	rec_rates->number_of_values=pop_size;
	gt_labels->number_of_values=pop_size;
	for (i=pop_size; i<number_of_individuals; i++)
	{
		for (w=0; w<number_of_words; w++) genotypes[index(i,w)]=NO_GENOTYPE;
	}
	//calculate its fitness and recombination rates
	calc_fit();
	calc_gt_labels();
	calc_rec();
	calc_stat();	//make sure everything is calculated
	if (HP_VERBOSE) cerr <<"done."<<endl;
	return 0;
}*/


/*
 * init with genotypes with 0
 */
int haploid_clone::init_genotypes(int n_o_genotypes)
{
	trait[0].coefficients_epistasis.clear();

	if (!mem)
	{
		cerr <<"haploid_clone::init_genotypes(): Allocate and set up first!\n";
		return HP_MEMERR;
	}
	if (n_o_genotypes>number_of_individuals)
	{
		cerr <<"haploid_clone::init_genotypes(): number of genotypes has to be smaller than max! got: "<<n_o_genotypes<<"!\n";
		return HP_BADARG;
	}
	pop_size=0;
	if (HP_VERBOSE) cerr <<"haploid_clone::init_genotypes(int n_o_genotypes, int no_copies) with population size "<<n_o_genotypes<<"...";
	generation=0;
	current_pop->clear();
	random_sample.clear();
	boost::dynamic_bitset<> tempgt(number_of_loci);
	tempgt.reset();			//set all bits to zero
	if (n_o_genotypes>=0){	//add n_o_genotypes copies of this genotype
		add_genotypes(tempgt,n_o_genotypes);
	}else{
		add_genotypes(tempgt,target_pop_size);
	}

	//calculate its fitness and recombination rates
	calc_fit();
	calc_stat();	//make sure everything is calculated
	if (HP_VERBOSE) cerr <<"done."<<endl;
	return 0;
}


/*
 * calculate the traits for all clones in the population
 * to be called when the fitness function is changed
 */
void haploid_clone::calc_everybodies_traits(){
	if (HP_VERBOSE) cerr<<"haploid_clone::calc_everybodies_traits()...";
	unsigned int ti;
	for (unsigned int i=0; i<current_pop->size(); i++)
	{
		calc_individual_fitness(&(*current_pop)[i]);
	}
}

/*
 * Calculate allele frequencies and mean and variance in fitness
 */
void haploid_clone::calc_stat()
{
	if (HP_VERBOSE) cerr<<"haploid_clone::calc_stat()...";
	calc_fitness_stat();
	int locus,cs;
	pop_size=0;
	for (locus=0; locus<number_of_loci; locus++) allele_frequencies[locus]=0;
	//loop over all clones
	for (unsigned int i=0; i<current_pop->size(); i++)
	{
		cs=(*current_pop)[i].clone_size;
		for (locus=0; locus<number_of_loci; locus++){	//add clone size to allele frequency of clone carries allele
			if ((*current_pop)[i].genotype[locus]) allele_frequencies[locus]+=cs;
		}
		pop_size+=cs;
	}
	//convert counts into frequencies
	for (locus=0; locus<number_of_loci; locus++)
	{
		allele_frequencies[locus]/=pop_size;
	}
	if (HP_VERBOSE) cerr<<"done.\n";
}

/*double haploid_clone::get_multi_point_frequency(vector <int> loci)
{
	if (HP_VERBOSE) cerr<<"haploid_clone::calc_multi_point_frequency()...";
	int *gt_mask=new int [number_of_words];
	int w,k,i;
	unsigned int  locus;
	double frequency=0;
	for (w=0; w<number_of_words; w++) gt_mask[w]=0;
	for (locus=0; locus<loci.size(); locus++)
	{
		w=locus_word(loci[locus]);
		k=locus_bit(loci[locus]);
		gt_mask[w]+=(1<<k);
		//cerr <<w<<" "<<k<<" "<<gt_mask[w]<<" ";
	}
	int temp;
	for (i=0; i<pop_size; i++)
	{
		temp=0;
		for (w=0; w<number_of_words; w++)
		{
			if ((genotypes[index(i,w)]&gt_mask[w])==gt_mask[w]) temp++;
		}
		if (temp==number_of_words) {
			frequency+=1.0;
		}
	}
	if (HP_VERBOSE) cerr<<"done.\n";
	delete [] gt_mask;
	//cerr <<" "<< frequency/pop_size<<endl;
	return frequency/pop_size;
}
*/

double haploid_clone::get_pair_frequency(int locus1, int locus2)
{
	if (HP_VERBOSE) cerr<<"haploid_clone::get_pair_frequency()...";
	unsigned int i;
	double frequency=0;
	for (i=0; i<current_pop->size(); i++)
		if ((*current_pop)[i].genotype[locus1] and (*current_pop)[i].genotype[locus2])
		        frequency += (*current_pop)[i].clone_size;
	if (HP_VERBOSE) cerr<<"done.\n";
	return frequency / pop_size;
}

vector <double>  haploid_clone::get_pair_frequencies(vector < vector <int> > *loci)
{
	if (HP_VERBOSE) cerr<<"haploid_clone::get_pair_frequencies()...";
	unsigned int i, pair;

	vector <double> freq;
	for (pair=0; pair<loci->size(); pair++) freq.push_back(0.0);
	for (i=0; i<current_pop->size(); i++)
	{
		for (pair=0; pair<loci->size(); pair++){
			if ((*current_pop)[i].genotype[(*loci)[pair][0]] and (*current_pop)[i].genotype[(*loci)[pair][1]]) freq[pair]+=(*current_pop)[i].clone_size;
		}
	}
	if (HP_VERBOSE) cerr<<"done.\n";
	for (pair=0; pair<loci->size(); pair++) freq[pair]/=pop_size;
	return freq;
}



/*
 * Do the evolution step
 */
int haploid_clone::evolve(){
	int err=0;
	random_sample.clear(); 		//discard the old random sample
	err+=select_gametes();		//select a new set of gametes (partitioned into sex and asex)
	err+=add_recombinants();	//do the recombination between pairs of sex gametes
	err+=swap_populations();	//make the new population the current population
	return err;
}

/*
 * for each clone, recalculate its fitness
 */
void haploid_clone::update_fitness(){
	for (unsigned int i=0; i<current_pop->size(); i++)
	{
		calc_fitness_from_traits(&(*current_pop)[i]);
	}
}

/*
 * Do the selection step
 */
int haploid_clone::select_gametes()
{
	//determine the current mean fitness, which includes a term to keep the population size constant
	double cpot=chemical_potential();
	//draw gametes according to parental fitness
	int os,o, nrec;
	pop_size=0;

	//to speed things up, reserve the expected amount of memory for sex gametes and the new population (+10%)
	sex_gametes.reserve(pop_size*outcrossing_probability*1.1);
	sex_gametes.clear();
	new_pop->reserve(current_pop->size()*1.1);
	new_pop->clear();
	for (unsigned int i=0; i<current_pop->size(); i++)
	{
		//poisson distributed random numbers -- mean exp(f)/bar{exp(f)}) (since death rate is one, the growth rate is (f-bar{f})
		if ((*current_pop)[i].clone_size>0){
			//the number of asex offspring of clone[i] is poisson distributed around e^{F-mF}(1-r)
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
	return 0;
}

/*
 * using the previously produced list of sex_gametes, pair them at random and mate
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

/*
 * After the new population is completely assembled, we have to make it the current
 * population. this is done simply by swaping pointers
 */
int haploid_clone::swap_populations(){
	//make the new generation the old one
	vector <gt> *temp_gt;
	temp_gt=current_pop;
	current_pop=new_pop;
	new_pop=temp_gt;

	generation++;
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
	gt offspring1(number_of_traits);
	gt offspring2(number_of_traits);
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

/*
 * choose a certain number of crossover points and produce a 0000111100011111101101 crossover pattern
 */
boost::dynamic_bitset<> haploid_clone::crossover_pattern()
{
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

/*
 * produce a random sample of genotypes (labeled by their clone of origin)
 * whcih is to be used for mutations and sampling.
 */
void haploid_clone::produce_random_sample(int size){
	if (HP_VERBOSE) cerr<<"haploid_clone::produce_random_sample(): size "<<size;

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

/*
 * return the clone index of a random genotype and remove it from the random_sample vector
 */
int haploid_clone::random_clone()
{
	int rclone;
	if (random_sample.size()>1){
		rclone=random_sample.back();
		random_sample.pop_back();
		return rclone;
	}else{	//if no genotypes left in sample, produce new sample
		produce_random_sample(min(pop_size, 1000));
		rclone=random_sample.back();
		random_sample.pop_back();
		return rclone;
	}
}


//mutate all loci
void haploid_clone::mutate()
{
	if (HP_VERBOSE)
		cerr <<"haploid_clone::mutate() .....";
	int i, actual_n_o_mutations,locus=0;
	produce_random_sample(number_of_loci*mutation_rate*pop_size*2);
	for (locus = 0; locus<number_of_loci; locus++) {
		//draw the number of mutation that are to happen at this locus
		actual_n_o_mutations = gsl_ran_poisson(evo_generator, mutation_rate*pop_size);
		//introduce these mutations one by one
		for (i = 0; i < actual_n_o_mutations; ++i) {
			flip_single_locus(random_clone(), locus);
		}
	}
	if (HP_VERBOSE)
		cerr <<"done";
}

//flip a spin at locus in random individual.
int haploid_clone::flip_single_locus(int locus)
{
	int clone = random_clone();
	flip_single_locus(clone, locus);
	return clone;
}


//flip a spin (locus) in individual, assign new fitness and recombination rate.
void haploid_clone::flip_single_locus(int clone, int locus)
{
	//produce new genotype
	gt tempgt(number_of_traits);
	tempgt.genotype.resize(number_of_loci);
	tempgt.genotype= (*current_pop)[clone].genotype;
	//new clone size == 1, old clone reduced by 1
	tempgt.clone_size=1;
	(*current_pop)[clone].clone_size--;
	//flip the locus in new clone and calculate fitness
	tempgt.genotype.flip(locus);
	calc_individual_fitness(&tempgt);
	//add clone to current population
	current_pop->push_back(tempgt);

	if (HP_VERBOSE) cerr <<"subpop::flip_single_spin(): mutated individual in clone "<<clone<<" at locus "<<locus<<endl;
}

/*
 * calculate the traits for the genotype and use them to calculate the fitness
 */
void haploid_clone::calc_individual_fitness(gt *tempgt)
{
	//calculate the new fitness value of the mutant
	for (int t=0; t<number_of_traits; t++){
		tempgt->trait[t] = trait[t].get_func(&(tempgt->genotype));
	}
	calc_fitness_from_traits(tempgt);
}

/*
 * convenience function, returns a genoytpe as string
 */
string haploid_clone::get_genotype_string(int i){
	string gt;
	boost::to_string((*current_pop)[i].genotype, gt);
	return gt;
}

/*
 * add the genotype specified by a bitset to the current population in in n copies
 */
void haploid_clone::add_genotypes(boost::dynamic_bitset<> genotype, int n)
{
	gt tempgt(number_of_traits);
	tempgt.genotype=genotype;
	tempgt.clone_size = n;
	calc_individual_fitness(&tempgt);
	current_pop->push_back(tempgt);
	pop_size+=n;
}

/*
 * return the chemical potential, i.e. the current mean fitness
 * plus term that causes the population size to relax to target_pop_size
 */
double haploid_clone::chemical_potential()
{
	calc_fitness_stat();
	return fitness_stat.mean+ (MIN(0.6931*(double(pop_size)/target_pop_size-1),2.0));
}

/*
 * go over every clone and calculate its fitness, followed by calculating the fitness statistics
 */
void haploid_clone::calc_fit()
{
	if (HP_VERBOSE)
		cerr <<"haploid_clone::calc_fit() .....";
	//recalculate the fitness of each clone
	update_fitness();
	calc_fitness_stat();
	if (HP_VERBOSE)
		cerr <<"done"<<endl;
}

/*
 * convenience function: return the fitness of the fittest clone
 */
double haploid_clone::get_max_fitness()
{
	double mf=(*current_pop)[0].fitness;
	for (unsigned int i=0; i<current_pop->size(); i++)
	{
		if (mf<(*current_pop)[i].fitness) mf=(*current_pop)[i].fitness;
	}
	return mf;
}


/*
 * calculate the statistics of fitness in the population
 */
void haploid_clone::calc_fitness_stat(){
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
		fitness_stat.mean+=temp*csize;
		fitness_stat.variance+=temp*temp*csize;
		pop_size+=csize;
	}
	if (pop_size==0){
		cerr <<"haploid_clone::calc_fitness_stat(): population extinct! clones: "<<current_pop->size()<<endl;
	}
	fitness_stat.mean/=pop_size;
	fitness_stat.variance/=pop_size;
	fitness_stat.variance-=fitness_stat.mean*fitness_stat.mean;
	if (HP_VERBOSE) {cerr <<"done"<<endl;}
}

/*
 * calculate trait statistics in the population and the trait covariances
 */
void haploid_clone::calc_trait_stat(){
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
int haploid_clone::print_allele_frequencies(ostream &out)
{
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
 * Calculate mean divergence from the [00...0] bitset
 */
double haploid_clone::divergence_mean()
{
	int n_sample = 1000;
	
	// Random number generator
	int seedtmp = time(NULL) + getpid() + generation;
	gsl_rng* int_gen = gsl_rng_alloc(RNG);
	gsl_rng_set(int_gen, seedtmp);

	// Cumulative partition the population according to clones
	vector <unsigned long> partition_cum = partition_cumulative();

	// Calculate divergence of the sample
	unsigned long NN = N(); 
	unsigned long c,tmp;
	boost::dynamic_bitset<> genotype;
	double divergence = 0.0;
	for (int i = 0; i < n_sample; i++) {
		c = gsl_rng_uniform_int(int_gen, NN);
		// What clones do they belong to?
		tmp = 0;
		while (c > partition_cum[tmp])	tmp++;
		c = tmp;
		// Calculate divergence
		genotype = (*current_pop)[c].genotype;
		for (int i = 0; i < genotype.size(); i++)
			divergence += genotype[i];
	}
	gsl_rng_free(int_gen);
	divergence /= n_sample;
	return divergence;
}

/*
 * Calculate diversity in the current population, i.e. Hamming distance between all pairs of sequences, and average.
 */
double haploid_clone::diversity_mean()
{
	int n_sample = 1000;
	
	// Random number generator
	int seedtmp = time(NULL) + getpid() + generation;
	gsl_rng* int_gen = gsl_rng_alloc(RNG);
	gsl_rng_set(int_gen, seedtmp);

	// Cumulative partition the population according to clones
	vector <unsigned long> partition_cum = partition_cumulative();

	// Calcolate random distances
	double diversity = 0.0;
	long tmp;
	unsigned long NN = N(); 
	unsigned long c, c1;
        c = c1 = 0;
	for (int i = 0; i < n_sample; i++) {
		while (c == c1) {
			c = gsl_rng_uniform_int(int_gen, NN);
			c1 = gsl_rng_uniform_int(int_gen, NN);
		}
		// What clones do they belong to?
		tmp = 0;
		while (c > partition_cum[tmp])	tmp++;
		c = tmp;
		tmp = 0;
		while (c1 > partition_cum[tmp])	tmp++;
		c1 = tmp;
		// Calculate distance if they belong to different clones
		if (c != c1 )
			diversity += distance_Hamming((*current_pop)[c].genotype,(*current_pop)[c1].genotype);
		c = c1 = 0;
	}
	gsl_rng_free(int_gen);
	diversity /= n_sample;
	return diversity;
}
/*
 * Calculate the hamming distance between two sequences (not normalized)
 */
unsigned int haploid_clone::distance_Hamming(boost::dynamic_bitset<> genotype, boost::dynamic_bitset<> genotype1)
{
	unsigned int d = 0;
	for (int i = 0; i != genotype.size(); i++)
		d += (genotype[i] != genotype1[i]);
	return d;
}

/*
 * Calculate the cumulative partition of sequences in the clones
 */
vector <unsigned long> haploid_clone::partition_cumulative()
{
	vector <unsigned long> partition_cum;
	unsigned long tmp;
	partition_cum.push_back((*current_pop)[0].clone_size);		
	for (int i = 1; i < get_number_of_clones(); i++) {
		partition_cum.push_back((*current_pop)[i].clone_size + partition_cum[i-1]);		
	}
	return partition_cum;
}

/*
 * function that reads the output of Hudson's ms and uses it to initialize the genotype distribution.
 * ms loci are fed into the genotype with a locus that is skipped. Each ms genotype is added multiplicity times
 */
int haploid_clone::read_ms_sample(istream &gts, int skip_locus, int multiplicity){
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
int haploid_clone::read_ms_sample_sparse(istream &gts, int skip_locus, int multiplicity, int distance){
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
				for (site=0; site<segsites and site*distance<L(); site++){
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
