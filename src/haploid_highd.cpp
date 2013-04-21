// vim: tabstop=8:softtabstop=8:shiftwidth=8:noexpandtab
/*
 * haploid_highd.cpp
 *
 *  Created on: Aug 28, 2008
 *      Author: neher
 *    Modified: Dec 13, 2011
 *   Committer: Fabio Zanini
 *
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
 */
#include <math.h>
#include "ffpopsim_highd.h"

/* Initialize the number of instances to zero */
size_t haploid_highd::number_of_instances = 0;

/**
 * @brief Default constructor
 *
 * The population objects are initialized, but neither the random generator is initialized nor the memory
 * for traits, gametes etc. is allocated. Please call set_up() before using the class.
 *
 * Note: The sequence is assumed to be linear (not circular). You can change this by hand if you wish so.
 */
haploid_highd::haploid_highd(int L_in, int rng_seed, int n_o_traits, bool all_polymorphic_in) {
	if (L_in < 1 or n_o_traits < 1) {
		if(HP_VERBOSE) cerr <<"haploid_highd::haploid_highd(): Bad Arguments! Both L and the number of traits must be larger or equal one."<<endl;
		throw HP_BADARG;
	}

	// Check that all_polymorphic is not used with more than one phenotypic trait
	if  (n_o_traits > 1 and all_polymorphic_in) {
		if(HP_VERBOSE) cerr <<"haploid_highd::haploid_highd(): Bad Arguments! Use multiple traits XOR all_polymorphic."<<endl;
		throw HP_BADARG;
	}

	// Check the bits per block (used in random reassortment)
	boost::dynamic_bitset<> temp;
	if (8*sizeof(long int) != temp.bits_per_block) {
		cerr <<"haploid_highd::haploid_highd(): haploid_highd requires sizeof(long int) to be equal to bits_per_block of boost::dynamic_bitset";
		throw HP_MEMERR;
	}

	// Set attributes
	number_of_loci = L_in;
	number_of_traits = n_o_traits;
	population_size = 0;
	number_of_clones = 0;
	mem = false;
	cumulants_mem = false;
	generation = -1;
	circular = false;
	carrying_capacity = 0;
	mutation_rate = 0;
	outcrossing_rate = 0;
	crossover_rate = 0;
	recombination_model = CROSSOVERS;
	fitness_max = HP_VERY_NEGATIVE;
	all_polymorphic=all_polymorphic_in;
	growth_rate = 2.0;

	//In case no seed is provided, get one from the OS
	seed = rng_seed ? rng_seed : get_random_seed();

	// Note: we should clean up the mess made by allocate_mem(). This requires more fine-grained
	// control than we currently have.
	int err = allocate_mem();
	if(err)	throw err;

	number_of_instances++;
}

/**
 * @brief Destructor
 *
 * Memory is released here.
 */
haploid_highd::~haploid_highd() {
	free_mem();
	number_of_instances--;
}

/**
 * @brief Get a random seed from /dev/urandom
 *
 * @returns non-deterministic, random seed
 */
int haploid_highd::get_random_seed() {
	int seedtmp;
	ifstream urandom("/dev/urandom", ios::binary);
	if(urandom.bad()) {
		cerr<<"/dev/urandom gives bad stream, falling back to time + getpid + number_of_instances"<<endl;
		seedtmp = time(NULL) + getpid() + number_of_instances;
	} else {
		urandom.read(reinterpret_cast<char*>(&seedtmp),sizeof(seedtmp));
		urandom.close();
	}
	return seedtmp;
}

/**
 * @brief Allocate all the necessary memory, initialze the RNG
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_highd::allocate_mem() {
	if (mem) {
		cerr <<"haploid_highd::allocate_mem(): memory already allocated, freeing and reallocating ...!\n";
		free_mem();
	}

	if (HP_VERBOSE) cerr <<"haploid_highd::allocate_mem(): loci: "<<number_of_loci<<endl;

	//Random number generator
	evo_generator = gsl_rng_alloc(RNG);
	gsl_rng_set(evo_generator, seed);
	if (HP_VERBOSE) cerr <<"haploid_highd() random number seed: "<<seed<<endl;
	//allocate all the memory
	genome = new int [number_of_loci+1];					// aux array holding range(0,number_of_loci) used to draw crossover points
	for (int i = 0; i < number_of_loci; i++) genome[i] = i;
	crossovers= new int [number_of_loci];					// aux array holding crossover points
	rec_pattern.resize(number_of_loci, 0);

	if (HP_VERBOSE) cerr <<"allele frequencies...";
	allele_frequencies = new double [number_of_loci];
	gamete_allele_frequencies = new double [number_of_loci];		//allele frequencies after selection

	trait = new hypercube_highd [number_of_traits];				//genotype trait function
	trait_stat = new stat_t [number_of_traits];				//structure holding trait statistics
	trait_covariance = new double* [number_of_traits];
	trait_weights = new double [number_of_traits];
	//initialize trait functions
	for (int t = 0; t < number_of_traits; t++){
		trait[t].set_up(number_of_loci, gsl_rng_uniform_int(evo_generator, 1<<20));
		trait_covariance[t] = new double [number_of_traits];
		trait_weights[t] = 0;
	}
	trait_weights[0] = 1;							// only the first trait is not zero

	mem = true;								//set memory flag to true
	if (HP_VERBOSE) cerr <<"done.\n";
	return 0;
}

/**
 * @brief Releases memory during class destruction
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_highd::free_mem() {
	if (!mem) {
		cerr <<"haploid_highd::free_mem(): No memory allocated!\n";
		return HP_BADARG;
	} else {
		delete [] allele_frequencies;
		delete [] gamete_allele_frequencies;
		delete [] trait;
		mem = false;
		return 0;
	}
}

/**
 * @brief: allocates memory for a sufficient number of clones
 *
 * @params: number of empty clones required
 */
int haploid_highd::provide_at_least(int n) {
	//calculate the number of clones that need to be newly allocated. Allow for some slack
	//to avoid calling this too often
	int needed_gts = n - available_clones.size() + 100 + 0.1 * population.size();

	//allocate at the necessary memory
	if (needed_gts > 50) {
		if (HP_VERBOSE) {cerr <<"haploid_highd::provide_at_least() requested: "<<n<<" providing: "<<needed_gts<<" total number of clones prev. allocated: "<<population.size()<< " number of clones available prev.: "<<available_clones.size()<<endl;}
		population.reserve(population.size()+needed_gts);

		//dummy clone used to push into the population
		clone_t tempgt(number_of_traits);
		tempgt.genotype.resize(number_of_loci,0);
		tempgt.clone_size=0;
		for (int ii = 0; ii < needed_gts; ii++) {
			available_clones.push_back(population.size());
			population.push_back(tempgt);
		}
		available_clones.reserve(population.size());
		sort(available_clones.begin(), available_clones.end(), std::greater<int>());
		if (track_genealogy){
			genealogy.extend_storage(population.size());
		}
		if (HP_VERBOSE) cerr <<" total number of clones new:: "<<population.size()<< " number of clones available: "<<available_clones.size()<<endl;
	}
	return 0;
}


/**
 * @brief Initialize population in linkage equilibrium
 *
 * @param freq target allele frequencies
 * @param N_in number of individuals to be created
 *
 * @returns zero if successful, error codes otherwise
 *
 * *Note*: the carrying capacity is set equal to N_in if it is still unset.
 */
int haploid_highd::set_allele_frequencies(double* freq, unsigned long N_in) {
	if (HP_VERBOSE) cerr <<"haploid_highd::set_allele_frequencies(double* freq, int N_in)...";
	if (N_in <= 0)	{
		cerr <<"The number of genotypes has to be positive!"<<endl;
		return HP_BADARG;
	}
	allele_frequencies_up_to_date=false;
	//reset the ancestral states
	ancestral_state.assign(L(), 0);
	polymorphism.assign(L(), poly_t());

	// set the carrying capacity if unset
	if(carrying_capacity < HP_NOTHING)
			carrying_capacity = N_in;

	// reset the current population
	population.clear();
	available_clones.clear();
	if (track_genealogy) {
		genealogy.reset_but_loci();
	}

	population_size = 0;
	number_of_clones = 0;
	last_clone = 0;
	provide_at_least(N_in);
        // set the allele frequencies
	boost::dynamic_bitset<> tempgt(number_of_loci);
	random_sample.clear();	//and the random sample
	if (HP_VERBOSE) cerr <<"add "<<N_in<<" genotypes of length "<<number_of_loci<<"..."<<endl;
	for (size_t i = 0; i < N_in; i++) {
		tempgt.reset();
		// set all loci for this genotype
		for(int locus = 0; locus < number_of_loci; locus++) {
			if (gsl_rng_uniform(evo_generator)<freq[locus])
				tempgt.set(locus);
		}
		add_genotype(tempgt,1);	//add genotype with multiplicity 1
		if (HP_VERBOSE >= 2) cerr <<i<<" ";
	}

	// Calculate all statistics to be sure
	generation++;
	calc_stat();

	//add the current generation to the genealogies and prune (i.e. remove parts that do not contribute the present)
	if (track_genealogy){genealogy.add_generation(fitness_max);}

	if (HP_VERBOSE) cerr <<"done."<<endl;
	return 0;
}



/**
 * @brief Initialize the population with genotype counts
 *
 * @param gt vector of genotype_value_pair with genotypes and sizes
 * @params vector that specifies the ancestral state of the sample. 
 *
 * @returns 0 if successful, nonzero otherwise.
 *
 * *Note*: the population size is set as the total sum of counts of all genotypes.
 * The carrying capacity is also set to the same number if it is still unset.
 */
int haploid_highd::set_genotypes(vector <genotype_value_pair_t> gt) {
  vector <int> tmp_ancestral_state(L(), 0); 
  return set_genotypes_and_ancestral_state(gt, tmp_ancestral_state);
}


/**
 * @brief Initialize the population with genotype counts
 *
 * @param gt vector of genotype_value_pair with genotypes and sizes
 * @params vector that specifies the ancestral state of the sample. 
 *
 * @returns 0 if successful, nonzero otherwise.
 *
 * *Note*: the population size is set as the total sum of counts of all genotypes.
 * The carrying capacity is also set to the same number if it is still unset.
 */
int haploid_highd::set_genotypes_and_ancestral_state(vector <genotype_value_pair_t> gt, vector <int>anc_state) {
	if (HP_VERBOSE) cerr <<"haploid_highd::set_genotypes_and_ancestral_state(vector <genotype_value_pair_t> gt)...";

	allele_frequencies_up_to_date = false;
	//reset the ancestral states
	ancestral_state.assign(L(), 0);
	polymorphism.assign(L(), poly_t());

	// Clear population
	population.clear();
	available_clones.clear();
	if (track_genealogy) {
		genealogy.reset_but_loci();
	}

	population_size = 0;
	random_sample.clear();

	// Initialize the clones and calculate the population size
	population_size = 0;
	number_of_clones = 0;
	last_clone = 0;
	provide_at_least(gt.size());

	if (anc_state.size()==L()){
	  for (size_t locus=0; locus<L(); locus++){
		ancestral_state[locus]=anc_state[locus];
	  }
	}else{
	  cerr <<"haploid_highd::set_genotypes_and_ancestral_state: length of ancestral state vector must equal number of loci"<<endl;
	  return HP_BADARG;
	}
	for(size_t i = 0; i < gt.size(); i++) {
		add_genotype(gt[i].genotype, gt[i].val);
		population_size += gt[i].val;
	}

	// set the carrying capacity if unset
	if(carrying_capacity < HP_NOTHING){carrying_capacity = population_size;}

	// Calculate all statistics to be sure
	generation++;
	calc_stat();

	//add the current generation to the genealogies and prune (i.e. remove parts that do not contribute the present)
	if (track_genealogy){genealogy.add_generation(fitness_max);}

	if (HP_VERBOSE) cerr <<"done."<<endl;
	return 0;
}

/**
 * @brief Initialize a wildtype population (00...0)
 *
 * @param N_in number of individuals
 *
 * @returns 0 if successful, nonzero otherwise.
 *
 * *Note*: the carrying capacity is set equal to N_in if it is still unset.
 */
int haploid_highd::set_wildtype(unsigned long N_in) {
	if (HP_VERBOSE) cerr <<"haploid_highd::set_wildtype(unsigned long N_in)...";

	if(N_in == 0) {
		if (HP_VERBOSE) cerr<<"the desired population size must be at least 1."<<endl;
		return HP_BADARG;
	}
	allele_frequencies_up_to_date = false;
	//reset the ancestral states
	ancestral_state.assign(L(), 0);
	polymorphism.assign(L(), poly_t());

	// Clear population
	population.clear();
	available_clones.clear();
	if (track_genealogy) {
		genealogy.reset_but_loci();
	}

	population_size = 0;
	number_of_clones = 0;
	last_clone = 0;
	random_sample.clear();
	provide_at_least(10);

	// Initialize the clones and calculate the population size
	boost::dynamic_bitset<> wildtype(number_of_loci);
	add_genotype(wildtype, N_in);

	// set the carrying capacity if unset
	if(carrying_capacity < HP_NOTHING)
			carrying_capacity = N_in;

	// Calculate all statistics to be sure
	generation++;
	calc_stat();
	// add the current generation to the genealogies and prune,
	// i.e. remove parts that do not contribute the present.
	if (track_genealogy){genealogy.add_generation(fitness_max);}
	if (HP_VERBOSE) cerr <<"done."<<endl;
	return 0;
}

/**

 * @brief Designates as set of loci to have their genealogy tracked
 *
 * @params locus to be tracked
 */
int haploid_highd::track_locus_genealogy(vector <int> loci) {
	//Note: you must track genealogies BEFORE the population is set
	if((generation != -1) or (get_number_of_clones() > 0)){
		cerr <<"haploid_highd::track_locus_genealogy: you must track genealogies BEFORE the population is set: generation "<<generation<<"\tnumber of clones"<<get_number_of_clones()<<endl;
		return HP_EXTINCTERR;
	}
	
	track_genealogy=true;
	if(HP_VERBOSE){cerr <<"haploid_highd::track_locus_genealogy(vector <int> loci)... number of loci="<<loci.size();}
	genealogy.reset();
	for (unsigned int i=0; i<loci.size(); i++){
		genealogy.track_locus(loci[i]);
	}
	genealogy.extend_storage(population.size());
	if (HP_VERBOSE){cerr<<"done\n";}
	return 0;
}


/**
 * @brief calculate and store allele frequencies
 *
 * Note: the allele frequencies are available in the allele_frequencies attribute.
 */
void haploid_highd::calc_allele_freqs() {
	// TODO: all iterations over L could be optimized with pointers
	if (HP_VERBOSE) cerr<<"haploid_highd::calc_allele_freqs()...";
	double cs;
	population_size = 0;
	participation_ratio = 0.0;
	for (int locus = 0; locus < number_of_loci; locus++)
		allele_frequencies[locus] = 0.0;

	//loop over all clones
	int i = 0;
	for(vector<clone_t>::iterator pop_iter = population.begin(); (pop_iter != population.end()) && (i < last_clone+1); pop_iter++, i++) {
		cs = pop_iter->clone_size;
		if (cs > 0) {
			for (int locus=0; locus<number_of_loci; locus++)	//add clone size to allele frequency of clone carries allele
				if (pop_iter->genotype[locus])
					allele_frequencies[locus] += cs;
			population_size += cs;
			participation_ratio += (cs * cs);
		}
	}
	//convert counts into frequencies
	participation_ratio /= population_size;
	participation_ratio /= population_size;
	for (int locus = 0; locus < number_of_loci; locus++)
		allele_frequencies[locus] /= population_size;
	if (HP_VERBOSE) cerr<<"done.\n";
	allele_frequencies_up_to_date = true;
}

/**
 * @brief Get the joint frequency of two alleles
 *
 * @param locus1 position of the first allele.
 * @param locus2 position of the second allele.
 *
 * @returns the joint frequency of the two alleles.
 */
double haploid_highd::get_pair_frequency(int locus1, int locus2) {
	if (HP_VERBOSE) cerr<<"haploid_highd::get_pair_frequency()...";

	double frequency = 0;
	int i = 0;
	for(vector<clone_t>::iterator pop_iter = population.begin(); (pop_iter != population.end()) && (i < last_clone+1); pop_iter++, i++)
		if ((pop_iter->clone_size > 0) and pop_iter->genotype[locus1] and pop_iter->genotype[locus2])
			frequency += pop_iter->clone_size;
	frequency /= population_size;

	if (HP_VERBOSE) cerr<<"done.\n";
	return frequency;
}

/**
 * @brief Get the joint frequency of two alleles, for a vector of allele pairs
 *
 * @param loci pointer to a vector of allele pairs. Each element of loci must be a vector of length 2.
 *
 * @returns vector of joint frequencies.
 */
vector <double>  haploid_highd::get_pair_frequencies(vector < vector <int> > *loci) {
	if (HP_VERBOSE) cerr<<"haploid_highd::get_pair_frequencies()...";

	unsigned int pair;

	vector <double> freq (loci->size(), 0.0);
	for (pair = 0; pair < loci->size(); pair++)
		freq[pair] = get_pair_frequency((*loci)[pair][0], (*loci)[pair][1]);

	if (HP_VERBOSE) cerr<<"done.\n";
	return freq;
}


/**
 * @brief Evolve for some generations under the specified conditions
 *
 * The order of steps performed is
 * 1. selection
 * 2. recombination
 * 3. mutation
 * but should not be important except for extremely high selection coefficients and/or very short times,
 * for which the discrete nature of the Fisher-Wright model becomes relevant.
 *
 * @param gen number of generations.
 *
 * @returns zero if successful, error code in the faulty step otherwise
 *
 * Note: if an error in encountered, evolution is stopped after the function that created the problem.
 * Typical errors include extinction or the opposite, offspring explosion.
 */
int haploid_highd::evolve(int gen) {
	if (HP_VERBOSE) cerr<<"haploid_highd::evolve(int gen)...";

	int err=0, g=0;
	allele_frequencies_up_to_date = false;
	// calculate an effective outcrossing rate to include the case of very rare crossover rates.
	// Since a recombination without crossovers is a waste of time, we scale down outcrossing probability
	// and scale up crossover rate so that at least one crossover is guaranteed to happen.
	if (recombination_model==CROSSOVERS)
		outcrossing_rate_effective = outcrossing_rate * (1 - exp(-number_of_loci * crossover_rate));
	else
		outcrossing_rate_effective = outcrossing_rate;

	// evolve cycle
	while((err == 0) && (g < gen)) {
		if (HP_VERBOSE) cerr<<"generation "<<generation<<endl;
		random_sample.clear();			//discard the old random sample
		if(err==0) err=select_gametes();	//select a new set of gametes (partitioned into sex and asex)
		else if(HP_VERBOSE) cerr<<"Error in select_gametes()"<<endl;
		sort(available_clones.begin(), available_clones.end(), std::greater<int>()); //sort clones in order to use the first ones again and again
		if(err==0) err=add_recombinants();	//do the recombination between pairs of sex gametes
		else if(HP_VERBOSE) cerr<<"Error in recombine()"<<endl;
		if(err==0) err=mutate();		//mutation step
		else if(HP_VERBOSE) cerr<<"Error in mutate()"<<endl;
		random_sample.clear();			//discard the old random sample
		g++;
		generation++;

		//add the current generation to the genealogies and prune (i.e. remove parts that do not contribute the present.
		if (track_genealogy) genealogy.add_generation(fitness_max);

	}
	if (HP_VERBOSE) {
		if(err==0) cerr<<"done."<<endl;
		else cerr<<"error "<<err<<"."<<endl;
	}
	return err;
}

/**
 * @brief Generate offspring according to fitness (selection) and segregate some for sexual mating
 *
 * @returns zero if successful, error codes otherwise
 *
 * Random Poisson offspring numbers are drawn from all parents, proportionally to their fitness.
 * A fraction r of all offspring are designated for sexual reproduction, while the remainder
 * (1-r) reproduces by exact duplication (mutations are introduced later).
 * The population size relaxes to a carrying capacity, i.e. selection is soft but the population
 * size is not exactly fixed.
 */
int haploid_highd::select_gametes() {
	// TODO: spot redundant recombination events
	if (HP_VERBOSE) cerr<<"haploid_highd::select_gametes()...";

	//determine the current mean fitness, which includes a term to keep the population size constant
	double relaxation = relaxation_value();
	allele_frequencies_up_to_date = false;

	//draw gametes according to parental fitness
	double delta_fitness;
	int os,o, nrec=0;
	int err = 0;
	population_size = 0;
	
	//sanity check
	if (outcrossing_rate_effective>1 or outcrossing_rate_effective<-HP_NOTHING) {
		cerr <<"haploid_highd::select_gametes(): outcrossing_rate needs to be <=1 and >=0, got: "<<outcrossing_rate_effective<<'\n';
		return HP_BADARG;
	}
	//to speed things up, reserve the expected amount of memory for sex gametes and the new population (+10%)
	sex_gametes.clear();
	clones_needed_for_recombination.clear();
	clones_needed_for_recombination.reserve(number_of_clones*outcrossing_rate_effective*1.1);
	sex_gametes.reserve(population_size*outcrossing_rate_effective*1.1);
	int new_last_clone = 0;
	number_of_clones = 0;
	fitness_max = HP_VERY_NEGATIVE;
	unsigned int clone_index = 0;
	for(vector<clone_t>::iterator pop_iter = population.begin(); (pop_iter != population.end()) && (clone_index < (unsigned int)(last_clone + 1)); pop_iter++, clone_index++){
		//poisson distributed random numbers -- mean exp(f)/bar{exp(f)})
		if (pop_iter->clone_size > 0) {
			//the number of asex offspring of clone[i] is poisson distributed around e^F / <e^F> * (1-r)
			delta_fitness = pop_iter->fitness - relaxation;
			//if (HP_VERBOSE >= 2) cerr<<i<<": relative fitness = "<<delta_fitness<<", Poisson intensity = "<<(pop_iter->clone_size*exp(delta_fitness)*(1-outcrossing_rate_effective))<<endl;
			//draw the number of sexual offspring, add them to the list of sex_gametes one by one
			if (outcrossing_rate_effective > 0){
				nrec = gsl_ran_poisson(evo_generator, pop_iter->clone_size * exp(delta_fitness) * outcrossing_rate_effective);
				for(o=0; o<nrec; o++) sex_gametes.push_back(clone_index);
			}

			os = gsl_ran_poisson(evo_generator, pop_iter->clone_size * exp(delta_fitness) * (1 - outcrossing_rate_effective));
			if (os > 0) {
				// clone[i] to new_pop with os as clone size
				pop_iter->clone_size = os;
				population_size += os;
				fitness_max = fmax(fitness_max, pop_iter->fitness);
				new_last_clone = clone_index;
				number_of_clones++;

			} else {
				pop_iter->clone_size = 0;
				if (nrec == 0)
					available_clones.push_back(clone_index);
				else
					clones_needed_for_recombination.push_back(clone_index);
			}
			if (track_genealogy){
				for (unsigned int locus=0; locus<genealogy.loci.size(); locus++){
					add_clone_to_genealogy(locus, clone_index, clone_index, 0,number_of_loci,os,(os>0));
				}
			}
		}
	}
	last_clone = new_last_clone;
	if(population_size+sex_gametes.size() < 1) {
		err = HP_EXTINCTERR;
		if (HP_VERBOSE) cerr<<"error "<<err<<". The population went extinct!"<<endl;
	}
	if ((HP_VERBOSE) and (err==0)) cerr<<"done."<<endl;
	return err;
}

/**
 * @brief Cause a bottleneck in the population size
 *
 * @param size_of_bottleneck number of individuals to leave alive (approximate)
 *
 * @returns zero if successful, error codes otherwise
 *
 * The bottleneck is performed as follows. Each clone is traded in for a smaller one with size given by
 * a Poisson random number around the new population size times the clone frequency. This function should
 * therefore almost conserve genotype frequencies, except for rare genotypes that are lost.
 *
 */
int haploid_highd::bottleneck(int size_of_bottleneck) {
// TODO: this function should accept a gsl random distribution as optional argument for choosing how sharp
// the bottleneck should be (i.e., how large fluctuations around the expected frequency may be). However,
// this requires function pointers or templates or lambda functions, and might be a nightmare to code.
	double ostmp;
	unsigned int os;
	int err = 0;
	unsigned int old_size = population_size;
	if (HP_VERBOSE) cerr<<"haploid_highd::bottleneck()...";
	allele_frequencies_up_to_date = false;

	population_size = 0;

	// resample each clone according to Poisson with a expected size reduced by bottleneck/N_old
	// eep track of the maximal fitness
	fitness_max = HP_VERY_NEGATIVE;
	unsigned int clone_index = 0;
	for(vector<clone_t>::iterator pop_iter = population.begin(); pop_iter != population.end() and (clone_index < (unsigned int)(last_clone + 1)); pop_iter++, clone_index++) {
		ostmp = pop_iter->clone_size * size_of_bottleneck / double(old_size);
		os = gsl_ran_poisson(evo_generator, ostmp);
		if(os > 0) {
			pop_iter->clone_size = os;
			population_size += os;
			check_individual_maximal_fitness(*pop_iter);
		} else {
			pop_iter->clone_size = 0;
			available_clones.push_back(clone_index);
		}
	}

	if(population_size < 1) err = HP_EXTINCTERR;

	if (HP_VERBOSE) {
		if(err==0) cerr<<"done."<<endl;
		else if(err == HP_EXTINCTERR) cerr<<" The population went extinct!"<<endl;
	}

	return err;
}


/**
 * @brief Mutate random clones at random loci (all loci)
 *
 * Note: for efficiency reasons, if the mutation rate is zero, an if cycle skips the body altogether.
 * The user can therefore call this function safely even if in non-mutating populations, without
 * loss of performance.
 *
 * Function is now rewritten to allow for multiple mutants.
 * TODO: Poisson conditional on at least one
 * FIXME!: all_polymorphic assumes that all additive effects are set, i.e. that locus equals the index of the
 * coefficient in the vector of additive effects
 */
int haploid_highd::mutate() {
	if (HP_VERBOSE)	cerr <<"haploid_highd::mutate() ..."<<endl;

	vector <int> mutations;
	int tmp_individual=0, nmut=0;
	size_t mutant;
	allele_frequencies_up_to_date = false;
	int actual_n_o_mutations,actual_n_o_mutants;
	if (mutation_rate > HP_NOTHING and not all_polymorphic) {
		//determine the number of individuals that are hit by at least one mutation
		actual_n_o_mutants = gsl_ran_poisson(evo_generator, (1.0-exp(-mutation_rate*number_of_loci))*population_size);
		produce_random_sample(min(actual_n_o_mutants, population_size));

		//make sure enough empty clones are available to accomodate the new mutants
		provide_at_least(min(actual_n_o_mutants, population_size));

		//loop over the mutant individuals and introduce the mutations
		for (int individual = 0; individual != actual_n_o_mutants; individual++) {
			//determine the target clone
			mutant = random_clone();
			//determine the number of mutation it suffers. this should be Poisson conditional on having at least one
			//in practice the solution is fine but it is somewhat inaccurate for multiple mutations
			actual_n_o_mutations = gsl_ran_poisson(evo_generator, number_of_loci * mutation_rate)+1;
			//introduce the mutations, not that flip_single_locus returns the number of new mutant, which is fed back into
			//flip_single_locus to introduce the next mutation
			for (int i = 0; i != actual_n_o_mutations; i++)
				mutant=flip_single_locus(mutant, gsl_rng_uniform_int(evo_generator,number_of_loci));
		}
	} else if(all_polymorphic) {
		if(HP_VERBOSE) cerr <<"haploid_highd::mutate(): keeping all loci polymorphic"<<endl;
		calc_allele_freqs(); //calculate the allele frequencies
		nmut=0;
		for (int locus=0; locus<L(); locus++){	//loop over all loci
			if (fabs(2*allele_frequencies[locus]-1)>1-HP_NOTHING){	//spot fixed loci
				if ((ancestral_state[locus]==0 and (2*allele_frequencies[locus]-1)<0) or
                    (ancestral_state[locus]==1 and (2*allele_frequencies[locus]-1)>0))                
                {	//if they are in the ancestral state
					tmp_individual = flip_single_locus(locus);		//introduce new allele
					polymorphism[locus].birth = get_generation();
					polymorphism[locus].fitness = population[tmp_individual].fitness-fitness_stat.mean;
					polymorphism[locus].fitness_variance = fitness_stat.variance;
					nmut++;
				}else{	//if locus is in derived state, flip coefficient of trait zero
					trait[0].set_additive_coefficient(-trait[0].get_additive_coefficient(locus),locus,locus);
					fixed_mutations.push_back(polymorphism[locus]);
					fixed_mutations.back().sweep_time = get_generation() -fixed_mutations.back().birth;
					tmp_individual=flip_single_locus(locus);
					ancestral_state[locus]= (ancestral_state[locus]==0)?1:0;
					polymorphism[locus].birth = get_generation();
					polymorphism[locus].effect = (2*ancestral_state[locus]-1)*trait[0].get_additive_coefficient(locus);
					polymorphism[locus].fitness = population[tmp_individual].fitness;
					polymorphism[locus].fitness_variance = fitness_stat.variance;
					nmut++;
				}
			}
		}
		number_of_mutations.push_back(nmut);
		calc_stat();

	} else if(HP_VERBOSE) cerr <<"haploid_highd::mutate(): mutation rate is zero."<<endl;

	if (HP_VERBOSE)	cerr <<"done."<<endl;;
	return 0;
}


/**
 * @brief Flip a spin at a specific locus in random individual
 *
 * @param locus position of the locus to flip
 *
 * @returns index of the new clone
 *
 * Note: this function calls flip_single_locus(unsigned int clonenum, int locus).
 */
unsigned int haploid_highd::flip_single_locus(int locus) {
	if (available_clones.size() == 0)
		provide_at_least(1);
	return flip_single_locus(random_clone(), locus);
}


/**
 * @brief Flip a spin (locus) in individual
 *
 * @param clonenum the individual whose locus is being flipped
 * @param locus position of the locus to flip
 *
 * This function creates a new clone and adds it to the population,
 * and assigns it a fitness.
 *
 * Note: This might produce duplicate clones since the mutant clone produced might
 * already exist. Duplicates can be merged by the member unique_clones()
 */
unsigned int haploid_highd::flip_single_locus(unsigned int clonenum, int locus) {
	// produce new genotype
	int new_clone = available_clones.back();
	available_clones.pop_back();
	allele_frequencies_up_to_date = false;

	//copy old genotype
	population[new_clone].genotype = population[clonenum].genotype;
	// new clone size == 1, old clone reduced by 1
	population[new_clone].clone_size = 1;
	population[clonenum].clone_size--;
	// flip the locus in new clone
	population[new_clone].genotype.flip(locus);
	// calculate traits and fitness
	vector<int> diff(1, locus);
	for (int t = 0; t < number_of_traits; t++){
		population[new_clone].trait[t] = population[clonenum].trait[t] + get_trait_difference(population[new_clone], population[clonenum], diff, t);
	}
	calc_individual_fitness_from_traits(population[new_clone]);
	check_individual_maximal_fitness(population[new_clone]);

	//update the last clones that is to be tracked
	last_clone = (new_clone<last_clone)?last_clone:new_clone;

	// add clone to current population
	if (population[clonenum].clone_size == 0)
		available_clones.push_back(clonenum);
	else
		number_of_clones++;

	if (track_genealogy) {
		for (unsigned int genlocus=0; genlocus<genealogy.loci.size(); genlocus++) {
			add_clone_to_genealogy(genlocus, new_clone,
					genealogy.newGenerations[genlocus][clonenum].parent_node.index,
					genealogy.newGenerations[genlocus][clonenum].crossover[0],
					genealogy.newGenerations[genlocus][clonenum].crossover[1], 1, 1);
			genealogy.newGenerations[genlocus][clonenum].clone_size--;
		}
	}


	if (HP_VERBOSE >= 2) cerr <<"subpop::flip_single_spin(): mutated individual in clone "<<clonenum<<" at locus "<<locus<<endl;
	return new_clone;
}

/**
 * @brief Pair and mate sexual gametes
 *
 * @returns zero if successful, nonzero otherwise
 *
 * Using the previously produced list of sex_gametes, pair them at random and mate
 */
int haploid_highd::add_recombinants() {
	//construct new generation
	int n_sex_gam = sex_gametes.size();
	int parent1, parent2, err;
	if (HP_VERBOSE) cerr <<"haploid_highd::add_recombinants(): add "<<n_sex_gam<<" recombinants!\n";

	if (n_sex_gam > 1) {
		//sexual offspring -- shuffle the set of gametes to ensure random mating
		gsl_ran_shuffle(evo_generator, &sex_gametes[0], n_sex_gam, sizeof(int));
		//make sure they are in even number
		if(n_sex_gam % 2) {sex_gametes.pop_back(); n_sex_gam--;}
		provide_at_least(n_sex_gam);
		//cout <<n_sex_gam<<'\t'<<population_size<<'\t'<<outcrossing_rate_effective<<endl;
		for(vector<int>::iterator iter = sex_gametes.begin(); iter != sex_gametes.end(); iter++) {
			parent1 = *iter;
			iter++;
			parent2 = *iter;
			//The recombination function stores two new genotypes
			//in new_genotypes[new_population_size] and [new_population_size+1]
			err = recombine(parent1, parent2);
			if(err) throw HP_RUNTIMEERR;
		}
	}
	for (vector<int>::iterator c = clones_needed_for_recombination.begin(); c != clones_needed_for_recombination.end(); c++)
		available_clones.push_back(*c);
	clones_needed_for_recombination.clear();
	return 0;
}


/**
 * @brief Recombine two genotypes parent1 and parent2 to produce two new genotypes
 *
 * @param parent1 first parent
 * @param parent2 second parent
 *
 * @returns zero if successful
 *
 * The new genotypes are stored in new_pop at positions ng and ng+1.
 *
 */
int haploid_highd::recombine(int parent1, int parent2) {
	if(HP_VERBOSE >= 2) cerr<<"haploid_highd::recombine(int parent1, int parent2)... parent 1: "<<parent1<<" parent 2: "<<parent2<<endl;

	allele_frequencies_up_to_date = false;

	//depending on the recombination model, produce a map that determines which offspring
	//inherites which part of the parental genomes
	if (recombination_model==FREE_RECOMBINATION)
		reassortment_pattern();
	else if (recombination_model==CROSSOVERS)
		crossover_pattern();
	//else {rec_pattern.resize(number_of_loci);}

	// produce two new genoytes
	int offspring_num1 = available_clones.back();
	available_clones.pop_back();
	number_of_clones++;
	int offspring_num2 = available_clones.back();
	available_clones.pop_back();
	number_of_clones++;

	if(HP_VERBOSE >= 2) cerr<<"offpring 1: "<<offspring_num1<<" offpring 2: "<<offspring_num2<<endl;

	// assign the genotypes by combining the relevant bits from both parents
	population[offspring_num1].genotype = (population[parent1].genotype&rec_pattern) | (population[parent2].genotype&(~rec_pattern));
	population[offspring_num2].genotype = (population[parent2].genotype&rec_pattern) | (population[parent1].genotype&(~rec_pattern));
	// clone size of new genoytpes is 1 each
	population[offspring_num1].clone_size = 1;
	population[offspring_num2].clone_size = 1;
	// calculate traits and fitness
	calc_individual_traits(population[offspring_num1]);
	calc_individual_traits(population[offspring_num2]);
	calc_individual_fitness_from_traits(population[offspring_num1]);
	calc_individual_fitness_from_traits(population[offspring_num2]);
	check_individual_maximal_fitness(population[offspring_num1]);
	check_individual_maximal_fitness(population[offspring_num2]);

	last_clone = (offspring_num1<last_clone)?last_clone:offspring_num1;
	last_clone = (offspring_num2<last_clone)?last_clone:offspring_num2;

	//Check what's going on
	if(HP_VERBOSE >= 3) {
		cerr<<rec_pattern<<endl;
		cerr<<population[parent1].genotype<<endl;
		cerr<<population[parent2].genotype<<endl;
		cerr<<population[offspring_num1].genotype<<endl;
		cerr<<population[offspring_num2].genotype<<endl<<endl;
	}

	population_size+=2;

	if (track_genealogy) {
		for (unsigned int genlocus=0; genlocus<genealogy.loci.size(); genlocus++){
			int locus=genealogy.loci[genlocus];
			int brleft=locus, brright=locus;
			bool state = rec_pattern[locus];
			while (rec_pattern[brleft]==state and brleft>0){brleft--;}
			while (rec_pattern[brright]==state and brright<number_of_loci){brright++;}
			brright--;
			if (state==1){
				add_clone_to_genealogy(genlocus, offspring_num1,parent1, brleft, brright, 1, 1);
				add_clone_to_genealogy(genlocus, offspring_num2,parent2, brleft, brright, 1, 1);
			}else{
				add_clone_to_genealogy(genlocus, offspring_num2,parent1, brleft, brright, 1, 1);
				add_clone_to_genealogy(genlocus, offspring_num1,parent2, brleft, brright, 1, 1);
			}
		}
	}

	if(HP_VERBOSE >= 2) cerr<<"done."<<endl;
	return 0;
}

void haploid_highd::add_clone_to_genealogy(int locusIndex, int dest, int parent, int left, int right, int cs, int n){
	if (HP_VERBOSE) {
		cerr <<"haploid_highd::add_clone_to_genealogy(): dest:  "<<dest<<" parent: "<<parent<<"  "<<genealogy.newGenerations[locusIndex].size()<<endl;
		tree_key_t temp;
		temp.age=generation-1;
		temp.index=parent;
		if (genealogy.trees[locusIndex].check_node(temp)){
			cerr <<"haploid_highd::add_clone_to_genealogy(): parent node ok"<<endl;
		}else{
			cerr <<"haploid_highd::add_clone_to_genealogy(): parent node DOES NOT EXIST!"<<endl;
		}
	}
	genealogy.newGenerations[locusIndex][dest].parent_node.index=parent;
	genealogy.newGenerations[locusIndex][dest].parent_node.age=generation-1;
	genealogy.newGenerations[locusIndex][dest].own_key.index=dest;
	genealogy.newGenerations[locusIndex][dest].own_key.age=generation;
	genealogy.newGenerations[locusIndex][dest].fitness= population[dest].fitness;
	genealogy.newGenerations[locusIndex][dest].number_of_offspring=n;
	genealogy.newGenerations[locusIndex][dest].clone_size=cs;
	genealogy.newGenerations[locusIndex][dest].crossover[0]=left;
	genealogy.newGenerations[locusIndex][dest].crossover[1]=right;
	if (HP_VERBOSE) {
		cerr <<"haploid_highd::add_clone_to_genealogy(): done"<<endl;
	}
}

/**
 * @brief For each clone, recalculate its traits
 */
void haploid_highd::update_traits() {
	int i=0;
	for(vector<clone_t>::iterator pop_iter = population.begin(); (pop_iter != population.end()) && (i<last_clone+1); pop_iter++, i++)
		if (pop_iter->clone_size>0)
			calc_individual_traits(*pop_iter);
}

/**
 * @brief For each clone, update fitness assuming traits are already up to date
 */
void haploid_highd::update_fitness() {
	if(population.size() > 0) {
		fitness_max = HP_VERY_NEGATIVE;
		unsigned int i = 0;
		for(vector<clone_t>::iterator pop_iter = population.begin(); (pop_iter != population.end()) and (i < (unsigned int)(last_clone + 1)); pop_iter++, i++)
			if (pop_iter->clone_size > 0) {
				calc_individual_fitness_from_traits(*pop_iter);
				check_individual_maximal_fitness(*pop_iter);
			}
	}
}

/**
 * @brief Calculate trait and fitness statistics and allele frequences
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
void haploid_highd::calc_stat() {
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
void haploid_highd::calc_individual_traits(clone_t &tempgt) {
	for (int t = 0; t < number_of_traits; t++)
		tempgt.trait[t] = trait[t].get_func(tempgt.genotype);
}

/**
 * @brief Calculate trait difference between two clones
 *
 * @param tempgt1 first clone
 * @param tempgt2 second clone
 * @param diffpos positions at which the genotypes differ
 * @param traitnum number of the trait to calculate
 *
 * @returns vector of differences
 *
 * The operation is, for each trait, tempgt1 - tempgt2.
 */
double haploid_highd::get_trait_difference(clone_t &tempgt1, clone_t &tempgt2, vector<int>& diffpos, int traitnum) {
	return trait[traitnum].get_func_diff(tempgt1.genotype, tempgt2.genotype, diffpos);
}

/**
 * @brief Calculate fitness from traits of the chosen clone
 *
 * @param tempgt clone whose fitness is to be calculated
 *
 * This function is linear in the traits with weights equal to trait_weights.
 * By default, only the first weight is different from zero.
 */
void haploid_highd::calc_individual_fitness_from_traits(clone_t &tempgt) {
	tempgt.fitness = trait_weights[0] * tempgt.trait[0];
	for (int t = 1; t < number_of_traits; t++)
		tempgt.fitness += trait_weights[t] * tempgt.trait[t];
}

/**
 * @brief Calculate fitness of a particular clone
 *
 * @param tempgt clone whose fitness is being calculated
 *
 * Note: this function also updates the traits information for the same clone, because the
 * phenotype is needed to calculate fitness. If you have already calculated the traits,
 * you can rely calc_individual_fitness_from_traits.
 */
void haploid_highd::calc_individual_fitness(clone_t &tempgt) {
	//calculate the new fitness value of the mutant
	calc_individual_traits(tempgt);
	calc_individual_fitness_from_traits(tempgt);
	//FIXME: why is this commented?
	//check_individual_maximal_fitness(tempgt);
}

/**
 * @brief Choose a number of crossover points and produce a crossover pattern
 *
 * @returns crossover pattern
 *
 * A typical crossover pattern would be 0000111100011111101101.
 */
void haploid_highd::crossover_pattern() {
	if (HP_VERBOSE) cerr<<"haploid_highd::crossover_pattern() "<<"...";

	int n_o_c = 0;
	vector <int> crossover_points;
	double total_rec = number_of_loci * crossover_rate;

	//TODO this should be poisson conditional on having at least one
	if (total_rec < 0.1) n_o_c=1;
	else while (n_o_c == 0) n_o_c = gsl_ran_poisson(evo_generator,total_rec);

	//for circular chromosomes make sure there is an even number of crossovers
	if (circular) {
		n_o_c *= 2;
		n_o_c = (n_o_c < number_of_loci)?n_o_c:number_of_loci;	//make sure there are fewer xovers than loci
		crossover_points.resize(n_o_c);
		//choose xovers at random from the genome label list
		//choose is expensive. could be replaced by simply random number followed by sorting
		gsl_ran_choose(evo_generator,(void*) &crossover_points[0],n_o_c,genome,number_of_loci,sizeof(int));
		for(vector<int>::iterator cp_iter = crossover_points.begin(); cp_iter != crossover_points.end(); cp_iter++)
			(*cp_iter)++; //increase all points by since crossover is after the selected locus
	} else {
		n_o_c = (n_o_c < number_of_loci)?n_o_c:(number_of_loci - 1);
		crossover_points.resize(n_o_c);
		for(vector<int>::iterator cp_iter = crossover_points.begin(); cp_iter != crossover_points.end(); cp_iter++)
			(*cp_iter) = gsl_rng_uniform_int(evo_generator,number_of_loci - 1) + 1;
		sort(crossover_points.begin(), crossover_points.end());
	}

	bool origin = true;
	rec_pattern.clear();
	//start with an empty bitset and extend to crossovers[c] with origing =0,1
	if (HP_VERBOSE>2) cerr<<" n_o_c: "<<n_o_c<<" origin "<<origin<<endl;
	for(vector<int>::iterator cp_iter = crossover_points.begin(); cp_iter != crossover_points.end(); cp_iter++) {
		if (HP_VERBOSE>2) {
			cerr<<" xo: "<<(*cp_iter)<<" origin "<<origin<<endl;
			cerr<<rec_pattern<<endl;
		}
		rec_pattern.resize((*cp_iter), origin);
		origin = !origin; //toggle origin
	}
	//extend to full length
	if (HP_VERBOSE>2) {
		cerr<<" xo: "<<" origin "<<origin;
		cerr<<rec_pattern<<endl;
	}
	rec_pattern.resize(number_of_loci, origin);

	if (HP_VERBOSE) {
		cerr<<rec_pattern<<endl;
		cerr <<"..done"<<endl;
	}
	return;
}

/**
 * @brief Produce a random reassortement pattern
 *
 * @returns reassortement pattern
 */
void haploid_highd::reassortment_pattern() {
	//the blocks of the bitset are to long for the rng, hence divide them by 4
	//TODO: this looks quite suspicious
	int bpblock = rec_pattern.bits_per_block / 4;
	long unsigned int temp_rec_pattern=0;
	int bits_left = number_of_loci, still_to_append;

	//set the random bitset with bpblock at a time
	rec_pattern.clear();
	while(bits_left >= 4 * bpblock) {
		temp_rec_pattern=gsl_rng_uniform_int(evo_generator, 1<<bpblock);
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<bpblock;
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<(2*bpblock);
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<(3*bpblock);
		//cout <<bpblock<<"  "<<temp_rec_pattern<<endl;
		rec_pattern.append(temp_rec_pattern);
		bits_left-=4*bpblock;
	}

	//set the remaining bits in blocks of bpblock and the remainder until nothing is left
	still_to_append = bits_left;
	if (bits_left >= 3 * bpblock) {
		temp_rec_pattern=gsl_rng_uniform_int(evo_generator, 1<<bpblock);
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<bpblock;
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<(2*bpblock);
		bits_left-=3*bpblock;
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bits_left)<<(3*bpblock);
	} else if (bits_left >= 2 * bpblock) {
		temp_rec_pattern=gsl_rng_uniform_int(evo_generator, 1<<bpblock);
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bpblock)<<bpblock;
		bits_left-=2*bpblock;
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bits_left)<<(2*bpblock);
	} else if (bits_left >= bpblock) {
		temp_rec_pattern=gsl_rng_uniform_int(evo_generator, 1<<bpblock);
		bits_left-=bpblock;
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bits_left)<<bpblock;
	} else
		temp_rec_pattern+=gsl_rng_uniform_int(evo_generator, 1<<bits_left);

	while (still_to_append) {
		still_to_append--;
		rec_pattern.push_back(((temp_rec_pattern&(1<<still_to_append))>0));
	}
}

/**
 * @brief Produce and store a random sample of the population for stochastic processes
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
void haploid_highd::produce_random_sample(int size) {
	if (HP_VERBOSE) cerr<<"haploid_highd::produce_random_sample(int): size "<<size<<"...";

	random_sample.clear();
	random_sample.reserve((size+50)*1.1);
	int thechosen, o;
	double frac = 1.1*(size+50)/population_size;
	//loop over all clones and choose a poisson distributed number of genoytpes
	unsigned int i = 0, cs;
	for(vector<clone_t>::iterator pop_iter = population.begin(); pop_iter != population.end() && (i < (unsigned int)(last_clone + 1)); pop_iter++, i++) {
		cs= pop_iter->clone_size;
		if (cs > 0) {
			thechosen = gsl_ran_poisson(evo_generator, frac*cs);
			//make sure it is not larger than the clone itself.
			thechosen = (pop_iter->clone_size < thechosen)?(cs):thechosen;
			//add each of the chosen individually to the random_sample vector
			if (thechosen) for (o = 0; o < thechosen; o++) random_sample.push_back(i);
		}
	}
	gsl_ran_shuffle(evo_generator, &random_sample[0], random_sample.size(), sizeof(int));
	if (HP_VERBOSE) cerr<<"done"<<endl;
}

/**
 * @brief Get a random clone from the population
 *
 * @returns the index of the random clone
 *
 * The probability density function from which the individual is chosen is flat over the
 * population (larger clones are proportionally more likely to be returned here).
 *
 * Note: for efficiency reasons, the class keeps a storage of random clone indices.
 * If you need much more than, say, 1000 random clones at a time, please call
 * produce_random_sample with the required size in advance.
 */
int haploid_highd::random_clone() {
	int rclone;
	int size = 1000;
	if (random_sample.size() > 1) {
		rclone=random_sample.back();
		random_sample.pop_back();
		return rclone;
	} else {
		//if no genotypes left in sample, produce new sample
		produce_random_sample(min(population_size, size));
		rclone=random_sample.back();
		random_sample.pop_back();
		return rclone;
	}
}

/**
 * @brief Sample random individuals from the population
 *
 * @param n_o_individuals number of individuals to sample
 * @param sample pointer to vector where to put the result
 *
 * @returns zero if successful, nonzero otherwise
 *
 * The results are not returned as a vector for performance reasons, as one might
 * want to get a lot of random clones. *sample may be not empty (but must be allocated).
 * In any case, clone numbers of the sampled individuals are appended to *sample. Hence,
 * you can use this function iteratively (although there might not be a good reason to
 * do so).
 */
int haploid_highd::random_clones(unsigned int n_o_individuals, vector <int> *sample) {
	sample->reserve(n_o_individuals);
	for(size_t i=0; i< n_o_individuals; i++)
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
void haploid_highd::add_genotype(boost::dynamic_bitset<> genotype, int n) {
	if(n > 0) {
		allele_frequencies_up_to_date = false;
		if (available_clones.size() == 0)
			provide_at_least(1);
		int new_gt = available_clones.back();
		available_clones.pop_back();

		population[new_gt].genotype = genotype;
		population[new_gt].clone_size = n;
		calc_individual_traits(population[new_gt]);
		calc_individual_fitness_from_traits(population[new_gt]);
		check_individual_maximal_fitness(population[new_gt]);

		population_size += n;
		last_clone = (new_gt < last_clone)?last_clone:new_gt;
		number_of_clones++;

		if (track_genealogy) {
			node_t leaf;
			leaf.fitness = population[new_gt].fitness;
			leaf.own_key.age=generation;
			leaf.own_key.index=new_gt;
			leaf.number_of_offspring = 1;
			leaf.clone_size = n;
			leaf.crossover[0]=0;
			leaf.crossover[1]=number_of_loci;
			for (unsigned int locusIndex=0; locusIndex<genealogy.loci.size(); locusIndex++){
				leaf.parent_node = genealogy.trees[locusIndex].get_MRCA();
				genealogy.newGenerations[locusIndex][new_gt]=leaf;
			}
		}
	}
}

/**
 * @brief Get the log of the exp-average fitness plus relaxation term
 *
 * @returns baseline relative to which growth rates are measured,
 * i.e. the current mean fitness plus the term that causes the population size to relax to carrying_capacity
 */
double haploid_highd::relaxation_value() {
	if (HP_VERBOSE) cerr <<"haploid_highd::relaxation_value()...";

	double logmean_expfitness = get_logmean_expfitness();
	// the second term is the growth rate when we start from N << carrying capacity
	double relax = logmean_expfitness + (fmin(log(growth_rate)*(double(population_size) / carrying_capacity - 1), 2.0)) + fitness_max;
	if (HP_VERBOSE)	cerr<<"log(<exp(F-Fmax)>) = "<<logmean_expfitness<<"... relaxation value = "<<relax<<"...done."<<endl;
	return relax;
}

/**
 * @brief Calculate and store fitness population statistics
 *
 * The fitness statistics are stored in fitness_stat.
 *
 * *Note*: this function assumes that fitness is up to date. If you are not sure, call update_traits() and update_fitness() first.
 */
void haploid_highd::calc_fitness_stat() {
	if (HP_VERBOSE) {cerr <<"haploid_highd::calc_fitness_stat()...";}

	double temp;
	int csize;
	fitness_stat.mean = 0;
	fitness_stat.variance = 0;
	population_size = 0;
	//loop over clones and add stuff up
	unsigned int i = 0;
	for(vector<clone_t>::iterator pop_iter = population.begin(); pop_iter != population.end() && (i < (unsigned int)(last_clone + 1)); pop_iter++, i++) {
		csize = pop_iter->clone_size;
		if (csize > 0) {
			temp = pop_iter->fitness;
			fitness_stat.mean += temp * csize;
			fitness_stat.variance += temp * temp * csize;
			population_size += csize;
		}
	}

	if (population_size == 0)
		cerr <<"haploid_highd::calc_fitness_stat(): population extinct! clones: "<<population.size()<<endl;

	//if (HP_VERBOSE) {cerr <<"pop size: "<<population_size<<", sum of fitnesses: "<<fitness_stat.mean<<"...";}
	fitness_stat.mean /= population_size;
	fitness_stat.variance /= population_size;
	fitness_stat.variance -= fitness_stat.mean * fitness_stat.mean;

	if (HP_VERBOSE) cerr <<"done."<<endl;
}


/**
 * @brief Get the population exp-average of fitness, used for keeping the population size fixed
 *
 * @returns the population exp-average of fitness
 *
 * Mathematically, this is \f$ \log\left( \left< e^(F-F_{max}) \right> \right) \f$. The baseline is \f$ F_{max} \f$
 * in order to avoid exponentiating large numbers.
 *
 */
double haploid_highd::get_logmean_expfitness() {
	if (HP_VERBOSE) cerr <<"haploid_highd::get_logmean_expfitness()...";

	double logmean_expfitness = 0;
	//loop over clones and add stuff up
	unsigned int i=0;
	for(vector<clone_t>::iterator pop_iter = population.begin(); pop_iter != population.end() && (i < (unsigned int)(last_clone + 1)); pop_iter++,i++)
		if (pop_iter->clone_size > 0)
			logmean_expfitness += pop_iter->clone_size * exp(pop_iter->fitness - fitness_max);
	logmean_expfitness /= population_size;
	logmean_expfitness = log(logmean_expfitness);
	if (HP_VERBOSE) cerr <<"done."<<endl;
	return logmean_expfitness;
}


/**
 * @brief Calculate and store trait population statistics and covariances
 *
 * The traits statistics are stored in traits_stat.
 *
 * *Note*: this function assumes that traits are up to date. If you are not sure, call update_traits() first.
 */
void haploid_highd::calc_trait_stat() {
	if (HP_VERBOSE) {cerr <<"haploid_highd::calc_trait_stat()...";}
	double temp,temp1;
	int t,t1, csize;
	population_size=0;

	// reset class attributes
	for(t = 0; t < number_of_traits; t++) {
		trait_stat[t].mean = 0;
		trait_stat[t].variance = 0;
		for(t1=0; t1<number_of_traits; t1++)
			trait_covariance[t][t1] = 0;
	}

	//loop over clones and add stuff up
	unsigned int i=0;
	for(vector<clone_t>::iterator pop_iter = population.begin(); pop_iter != population.end() && (i < (unsigned int)(last_clone + 1)); pop_iter++,i++) {
		csize = pop_iter->clone_size;
		if (csize>0) {
			for(t = 0; t < number_of_traits; t++) {
				temp = pop_iter->trait[t];
				trait_stat[t].mean += temp * csize;
				trait_stat[t].variance += temp * temp * csize;
				for(t1 = 0; t1 < number_of_traits; t1++) {
					temp1 = pop_iter->trait[t1];
					trait_covariance[t][t1] += temp * temp1 * csize;
				}
			}
			population_size += csize;
		}
	}
	//complain if population went extinct
	if (population_size == 0)
		cerr <<"haploid_highd::calc_trait_stat(): population extinct! clones: "<<population.size()<<endl;

	//normalize the means and variances
	for(t = 0; t < number_of_traits; t++) {
		trait_stat[t].mean /= population_size;
		trait_stat[t].variance /= population_size;
		trait_stat[t].variance -= trait_stat[t].mean * trait_stat[t].mean;
	}
	//calculate the covariances
	for(t = 0; t < number_of_traits; t++)
		for(t1 = 0; t1 < number_of_traits; t1++) {
			trait_covariance[t][t1] /= population_size;
			trait_covariance[t][t1] -= trait_stat[t].mean * trait_stat[t1].mean;
		}

	if (HP_VERBOSE) {cerr <<"done"<<endl;}
}

/**
 * @brief Print all allele frequencies into a stream provided
 *
 * @param out stream to put the allele frequencies (usually a file or stdout)
 *
 * @returns zero if successful, nonzero otherwise
 */
int haploid_highd::print_allele_frequencies(ostream &out) {
	if (out.bad()) {
		cerr <<"haploid_highd::print_allele_frequencies: bad stream\n";
		return HP_BADARG;
	}

	calc_stat();
	out <<setw(10)<<generation;
	for (int l=0; l<number_of_loci; l++)
		out<<setw(15)<<allele_frequencies[l];
	out <<endl;

	return 0;
}

/**
 * @brief Read the output of Hudson's ms and use it to initialize the genotype distribution
 *
 * @param gts genotypes output of _ms_
 * @param skip_locus position of the locus to be skipped
 * @param multiplicity number of times each genotype is added
 *
 * @returns zero if successful, error codes otherwise
 *
 * ms loci are fed into the genotype with a locus that is skipped. Each ms genotype is added multiple times.
 *
 */
int haploid_highd::read_ms_sample(istream &gts, int skip_locus, int multiplicity) {
	if (gts.bad()) {
		cerr<<"haploid_highd::read_ms_sample(): bad stream!\n";
		return HP_BADARG;
	}

	allele_frequencies_up_to_date = false;
	//line buffer to read in the ms input
	char *line = new char [2*number_of_loci+5000];
	bool found_gt = false;
	string header;
	int count = 0;
	int segsites, site, locus;
	segsites = 0;

	//new genotype to be read in from ms
	boost::dynamic_bitset<> newgt(number_of_loci);
	//reset population
	population.clear();
	random_sample.clear();
	population_size = 0;
	if (mem) {
		set_wildtype(0);
		population_size = 0;
		//loop over each line of the ms out-put file
		while(gts.eof() == false) {
			gts.get(line, 2*number_of_loci+5000);
			gts.get();
			//skip over empty lines
			while (gts.peek() == '\n')
				gts.get();
			//cout <<count<<"  "<<gt<<" "<<line<<endl;
			count++;
			//if the first genotype has been found
			if (found_gt and line[0] != '\0') {
				newgt.reset();
				//go over the line and assing loci, skip over the "skip_locus"
				for (site = 0; site < segsites; site++)
					if (line[site] == '1') {
						if (site < skip_locus) {
							locus=site;
							newgt.set(locus);
						} else {
							locus=site+1;
							newgt.set(locus);
						}
					}
				add_genotype(newgt, multiplicity);

			} else { //reading the header
				header.assign(line);
				if (header.compare(0,2, "//") == 0) {
					gts.get(line, 2*number_of_loci);
					gts.get();
					while (gts.peek() == '\n')
						gts.get();
					cerr <<count<<"  "<<found_gt<<" "<<line<<endl;
					header.assign(line);
					segsites = atoi(header.substr(9,header.size()-9).c_str());
					gts.get(line, 2*number_of_loci);
					gts.get();
					while (gts.peek() == '\n')
						gts.get();
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
 * @brief Read the output of Hudson's ms and use it to initialize the genotype distribution
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
int haploid_highd::read_ms_sample_sparse(istream &gts, int skip_locus, int multiplicity, int distance) {
	if (gts.bad()) {
		cerr<<"haploid_highd::read_ms_sample(): bad stream!\n";
		return HP_BADARG;
	}
	//line buffer to read in the ms input
	char *line = new char [2*number_of_loci+5000];
	bool found_gt = false;
	string header;
	int count = 0;
	int segsites, site, locus;
	segsites = 0;
	allele_frequencies_up_to_date = false;

	//new genotype to be read in from ms
	boost::dynamic_bitset<> newgt(number_of_loci);
	//reset population
	population.clear();
	random_sample.clear();
	population_size = 0;
	if (mem) {
		set_wildtype(0);
		population_size=0;
		//loop over each line of the ms out-put file
		while(gts.eof() == false) {
			gts.get(line, 2*number_of_loci+5000);
			gts.get();
			//skip over empty lines
			while (gts.peek() == '\n')
				gts.get();

			//cout <<count<<"  "<<gt<<" "<<line<<endl;
			count++;
			//if the first genotype has been found
			if (found_gt and line[0] != '\0') {
				newgt.reset();
				//go over the line and assing loci, skip over the "skip_locus"
				for (site = 0; site < segsites and site * distance < get_number_of_loci(); site++)
					if (line[site] == '1') {
					    locus = site * distance;
					    if (locus != skip_locus)
						    newgt.set(locus);
					}
				add_genotype(newgt, multiplicity);
			} else { //reading the header
				header.assign(line);
				if (header.compare(0,2, "//") == 0) {
					gts.get(line, 2*number_of_loci);
					gts.get();
					while (gts.peek() == '\n')
						gts.get();
					cerr <<count<<"  "<<found_gt<<" "<<line<<endl;
					header.assign(line);
					segsites = atoi(header.substr(9,header.size()-9).c_str());
					gts.get(line, 2*number_of_loci);
					gts.get();
					while (gts.peek() == '\n')
						gts.get();
					cerr <<count<<"  "<<found_gt<<" "<<line<<endl;
					found_gt = true;
				}
			}
		}
	}
	delete [] line;
	return 0;
}


/**
 * @brief Calculate Hamming distance between two sequences
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
int haploid_highd::distance_Hamming(boost::dynamic_bitset<> gt1, boost::dynamic_bitset<> gt2, vector <unsigned int *> *chunks, unsigned int every) {
	// check whether we have chunks at all
	if((!chunks) or (chunks->size() == 0)) {
		if(every!=1) return HP_BADARG;
		else return (gt1 ^ gt2).count();
	}

	unsigned int d = 0;
	unsigned int pos;

	// check that the chunks make sense
	if((every < 1) or (every >= (unsigned int)number_of_loci)) return HP_BADARG;
	for(vector<unsigned int *>::iterator ck_iter = chunks->begin(); ck_iter != chunks->end(); ck_iter++) {
		if((*ck_iter)[1] < (*ck_iter)[0]) return HP_BADARG;
		if((*ck_iter)[1] >= (unsigned int)number_of_loci) return HP_BADARG;
		if((*ck_iter)[0] >= (unsigned int)number_of_loci) return HP_BADARG;
		
		for(pos = (*ck_iter)[0]; pos < (*ck_iter)[1]; pos += every)
			d += (unsigned int)(gt1[pos] != gt2[pos]);
	}
	return d;
}


/**
 * @brief Calculate the cumulative partition of sequences into clones
 *
 * @param partition_cum the vector to be filled
 *
 * @returns vector of cumulative clone sizes
 *
 * *Example*: if there are three clones of sizes (100, 22, 3) this function will
 * fill the vector (100, 122, 125). The last element is of course the population
 * size See also get_population_size.
 *
 * *Note*: the vector is taken in input by reference for performance reasons, since it can get huge.
 */
int haploid_highd::partition_cumulative(vector <unsigned int> &partition_cum) {	
	if(population.size() == 0)
		return HP_EXTINCTERR;

	partition_cum.clear();
	partition_cum.push_back(population[0].clone_size);
	if(population.size() > 1) 
		for(vector<clone_t>::iterator pop_iter = population.begin() + 1; pop_iter != population.end(); pop_iter++)
			partition_cum.push_back(pop_iter->clone_size + partition_cum.back());
	return 0;
}

/**
 * @brief Calculate mean and variance of the divergence from the [00...0] bitset
 *
 * @param n_sample size of the statistical sample to use (the whole pop is often too large)
 *
 * @returns mean and variance of the divergence in a stat_t 
 */
stat_t haploid_highd::get_divergence_statistics(unsigned int n_sample) {
	stat_t div;
	unsigned int tmp;
	vector <int> clones;
	produce_random_sample(n_sample);
	random_clones(n_sample, &clones);

	for (size_t i = 0; i < n_sample; i++) {
		tmp = (population[clones[i]].genotype).count();
		div.mean += tmp;
		div.variance += tmp * tmp;
	}
	div.mean /= n_sample;
	div.variance /= n_sample;
	div.variance -= div.mean * div.mean;
	return div;
}

/**
 * @brief Calculate diversity in the current population (Hamming distance between pairs of sequences)
 *
 * @param n_sample size of the statistical sample to use (the whole pop is often too large)
 *
 * @returns mean and variance of the diversity in a stat_t
 */
stat_t haploid_highd::get_diversity_statistics(unsigned int n_sample) {
	stat_t div;
	unsigned int tmp;
	vector <int> clones1;
	vector <int> clones2;
	produce_random_sample(n_sample * 2);
	random_clones(n_sample, &clones1);
	random_clones(n_sample, &clones2);

	for (size_t i = 0; i < n_sample; i++) {
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
 * @brief Calculate histogram of fitness from traits
 *
 * @param hist pointer to the gsl_histogram to fill
 * @param bins number of bins in the histogram 
 *
 * *Note*: the output histogram might have less bins than requested if the sample size is too small.
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
int haploid_highd::get_fitness_histogram(gsl_histogram **hist, unsigned int bins, unsigned int n_sample) {
	if (HP_VERBOSE) cerr <<"haploid_highd::get_fitness_histogram()...";

	// Calculate fitness of the sample
	double fitnesses[n_sample];
	vector <int> clones;
	produce_random_sample(n_sample);
	random_clones(n_sample, &clones);
	for(size_t i = 0; i < n_sample; i++)
		fitnesses[i] = population[clones[i]].fitness;

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
 * @brief Get histogram of divergence from the [00...0] bitset
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
int haploid_highd::get_divergence_histogram(gsl_histogram **hist, unsigned int bins, vector <unsigned int *> *chunks, unsigned int every, unsigned int n_sample) {
	if (HP_VERBOSE) cerr <<"haploid_highd::get_divergence_histogram(gsl_histogram **hist, unsigned int bins, vector <unsigned int *> *chunks, unsigned int every, unsigned int n_sample)...";

	boost::dynamic_bitset<> gt_wt(number_of_loci);	// the [00...0] bitset
	int temp;
	unsigned int divs[n_sample];
	vector <int> clones;
	produce_random_sample(n_sample);
	random_clones(n_sample, &clones);
	for(size_t i = 0; i < n_sample; i++) {
		temp = distance_Hamming(gt_wt, population[clones[i]].genotype, chunks, every);
		// negative distances are error codes
		if(temp < 0) return temp;
		else divs[i] = temp;
	}

	// Prepare the histogram
	unsigned long dmax = *max_element(divs, divs + n_sample);
	unsigned long dmin = *min_element(divs, divs + n_sample);

	// Antialiasing
	unsigned int width, binsnew;
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
 * @brief Get histogram of diversity in the population (mutual Hamming distance)
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
int haploid_highd::get_diversity_histogram(gsl_histogram **hist, unsigned int bins, vector <unsigned int *> *chunks, unsigned int every, unsigned int n_sample) {
	if (HP_VERBOSE) {cerr <<"haploid_highd::get_diversity_histogram(gsl_histogram **hist, unsigned int bins, vector <unsigned int *> *chunks, unsigned int every, unsigned int n_sample)...";}
	int temp;
	unsigned int divs[n_sample];
	vector <int> clones1;
	vector <int> clones2;
	produce_random_sample(n_sample * 2);
	random_clones(n_sample, &clones1);
	random_clones(n_sample, &clones2);
	for(size_t i = 0; i < n_sample; i++) {
		temp = distance_Hamming(clones1[i], clones2[i], chunks, every);
		// negative distances are error codes
		if(temp < 0) return temp;
		else divs[i] = temp;
	}

	// Prepare the histogram
	unsigned long dmax = *max_element(divs, divs + n_sample);
	unsigned long dmin = *min_element(divs, divs + n_sample);

	// Aliasing
	unsigned int width, binsnew;
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
	gsl_histogram_scale(*hist, 1 / (double)n_sample);
	
	if (HP_VERBOSE) cerr<<"done.";
	return 0;
}

/**
 * @brief Remove duplicate clones.
 *
 * The library does not usually check whether two clones with the same genotype are present, but
 * this can happen in case of recurrent mutation. This function merges duplicates.
 *
 * *Note*: this is only needed for studying the clone structure. Evolution itself does not need to
 * make sure that clones are unique.
 */
void haploid_highd::unique_clones() {
	random_sample.clear();
	number_of_clones = 0;
	population_size = 0;
	int new_last_clone = 0;
	if(population.size() > 1) {
		// sort them O(nlog(n))
		sort(population.begin(), population.begin() + last_clone + 1);
		sort(available_clones.begin(), available_clones.end(), std::greater<int>());
		reverse(population.begin(), population.begin() + last_clone + 1);
		//available_clones.clear();
		while (available_clones.back() < last_clone+1) available_clones.pop_back();
		// merge clones with the same fitness and genotype
		vector<clone_t>::iterator pop_iter = population.begin();
		unsigned int i = 0;
		while (pop_iter->clone_size == 0) {
			pop_iter++;
			i++;
		}
		clone_t  *current_last_clone = &(*pop_iter);
		number_of_clones++;
		new_last_clone = i;
		population_size += pop_iter->clone_size;
		pop_iter++;
		i++;
		for(; pop_iter != population.end() && (i < (unsigned int)(last_clone + 1)); pop_iter++,i++) {
		//for(; pop_iter != population.end(); pop_iter++,i++) {
			if (pop_iter->clone_size > 0){
				if((*pop_iter) == (*current_last_clone)) {
					current_last_clone->clone_size += pop_iter->clone_size;
					available_clones.push_back(i);
					population_size += pop_iter->clone_size;
					pop_iter->clone_size=0;
				} else {
					current_last_clone= &(*pop_iter);
					number_of_clones++;
					population_size += pop_iter->clone_size;
					new_last_clone = i;
				}
			} else
				available_clones.push_back(i);
		}
		sort(available_clones.begin(), available_clones.end(), std::greater<int>());
		last_clone=new_last_clone;
	}
}

/**
 * @brief Obtain a list of the good clones.
 *
 * @returns vector of clone indices
 */
vector <int> haploid_highd::get_nonempty_clones() {
	vector <int> good;
	good.reserve(population.size());
	unsigned int i = 0;
	vector<clone_t>::iterator pop_iter = population.begin();
	for(; pop_iter != population.end() && (i < (unsigned int)(last_clone + 1)); pop_iter++,i++)
		if(pop_iter->clone_size)
			good.push_back(i);
	return good;
}

