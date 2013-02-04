/*
 * haploid_lowd.cpp
 *
 *  Created on: Jan 27, 2010
 *  Author: Richard Neher
 *  Modified by: Fabio Zanini
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
#include "ffpopsim_lowd.h"

/* Initialize the number of instances to zero */
size_t haploid_lowd::number_of_instances=0;

/**
 * @brief Default constructor
 *
 * @param L_in number of loci (at least 1)
 * @param rngseed seed for the random number generator. If this is zero, time(NULL)+getpid() is used.
 */
haploid_lowd::haploid_lowd(int L_in, int rng_seed) {
	if (L_in <1) {
		cerr <<"haploid_lowd::haploid_lowd(): Bad Arguments! L must be larger or equal one."<<endl;
		throw HG_BADARG;
	}

	// Set attributes
	number_of_loci = L_in;
	population_size = carrying_capacity = 0;
	recombination_model = FREE_RECOMBINATION;
	mem = false;
	outcrossing_rate = 1.0;
	circular = false;
	generation = 0;
	long_time_generation = 0;

	//In case no seed is provided, get one from the OS
	seed = rng_seed ? rng_seed : get_random_seed();

	// Note: we should clean up the mess made by allocate_mem(). This requires more fine-grained
	// control than we currently have.
	int err = allocate_mem();
	if(err)	throw err;
	
	number_of_instances++;
}

/**
 * @brief Default destructor
 *
 * Release memory.
 */
haploid_lowd::~haploid_lowd() {
	free_mem();
	number_of_instances--;
}

/**
 * @brief Get a random seed from /dev/urandom
 *
 * @returns non-deterministic, random seed
 */
int haploid_lowd::get_random_seed() {
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
 * @brief Allocate all the necessary memory, initialze the RNG.
 *
 * @returns zero if successful, error codes otherwise
 *
 * Set up the different hypercubes needed to store the fitness, population recombinants, and mutants.
 */
int haploid_lowd::allocate_mem() {
	if (mem) {
		if(HG_VERBOSE) cerr <<"haploid_lowd::allocate_mem(): Memory allocated already!"<<endl;
		return HG_BADARG;
	}

	int err = 0;
	rng = gsl_rng_alloc(RNG);
	gsl_rng_set(rng, seed);
	if (HG_VERBOSE) cerr <<"haploid_lowd() random number seed: "<<seed<<endl;
	err += fitness.set_up(number_of_loci, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
	err += population.set_up(number_of_loci, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
	err += mutants.set_up(number_of_loci, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
	err += recombinants.set_up(number_of_loci, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));

	mutation_rates = new double*[2]; //allocate backward and forward mutation rates arrays
	mutation_rates[0] = new double [number_of_loci]; //allocate forward mutation rates
	mutation_rates[1] = new double [number_of_loci]; //allocate backward mutation rates

	//set mutation rates to zero
	for (int fb = 0; fb < 2; fb++)
		for (int locus = 0; locus < number_of_loci; locus++)
			mutation_rates[fb][locus] = 0;
	if (err==0) {
		mem=true;
		return 0;
	}
	else
		return HG_MEMERR;
}


/**
 * @brief Releases memory during class destruction.
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_lowd::free_mem() {
	if (!mem) {
		if(HG_VERBOSE) cerr <<"haploid_lowd::free_mem(): No memory allocated!"<<endl;
		return HG_BADARG;
	}

	gsl_rng_free(rng);
	int err = free_recombination_mem();
	delete [] mutation_rates[1];
	delete [] mutation_rates[0];
	delete [] mutation_rates;
	mem=false;
	return err;
}

/**
 * @brief Allocate memory for recombination patterns
 *
 * @param rec_model recombination model to use (CROSSOVERS or SINGLE_CROSSOVER)
 *
 * @returns zero if successful, number of failed allocations otherwise
 *
 * *Note*: this function also sets the recombination_model flag.
 *
 * *Note*: in the SINGLE_CROSSOVER model, only recombination_patterns[0] exists,
 * and it stores the crossover rates.
 */
int haploid_lowd::allocate_recombination_mem(int rec_model) {
	if (HG_VERBOSE) cerr<<"haploid_lowd::allocate_recombination_mem()...";
	int err = 0;
	int spin, temp;
	int *nspins;	//temporary variables the track the number of ones in the binary representation of i
	if (rec_model == CROSSOVERS) {
		nspins=new int [1<<number_of_loci];
		if (nspins==NULL) {
			cerr<<"haploid_lowd::set_recombination_rates(): Can not allocate memory!"<<endl;
			return HG_MEMERR;
		}
		spin=-1;
		nspins[0]=0;
		//allocate space for all possible subsets of loci
		recombination_patterns = new double* [1<<number_of_loci];
		recombination_patterns[0] = new double [1];
		//loop over all possible locus subsets and allocate space for all
		//possible ways to assign the subset to father and mother (2^nspins)
		for (int i = 1; i < (1<<number_of_loci); i++) {
			// coefficients are sorted in a specific order
			// the order of coefficient k is 1+(the order of coefficient[k-2^spin])
			if (i==(1<<(spin+1))) spin++;
			temp = 1+nspins[i-(1<<spin)];
			nspins[i] = temp;

			//all possible ways to assign the subset to father and mother
			// for CROSSOVERS: 2^temp
			//     e.g. 000, 001, 010, 011, 100, 101, 110, 111
			recombination_patterns[i] = new double [(1<<temp)];
			if (recombination_patterns[i]==NULL) err += 1;
		}
		recombination_model = CROSSOVERS;
		delete [] nspins;

	} else if (rec_model == SINGLE_CROSSOVER) {
		recombination_patterns = new double* [1];
		recombination_patterns[0] = new double [number_of_loci];
		recombination_model = SINGLE_CROSSOVER;
	}

	if (HG_VERBOSE) cerr<<"...done."<<endl;
	return err;
}

/**
 * @brief Release memory for recombination patterns
 *
 * @returns zero if successful, error codes otherwise
 *
 * *Note*: this function is called every time the recombination rates are set with another
 * model for recombination than the previous one, and upon class destruction.
 */
int haploid_lowd::free_recombination_mem() {
	if (HG_VERBOSE) cerr<<"haploid_lowd::free_recombination_mem()...";

	if (recombination_model == CROSSOVERS) {
		for (int i = 0; i < (1<<number_of_loci); i++)
			delete [] recombination_patterns[i];
		delete [] recombination_patterns;

	} else if (recombination_model==SINGLE_CROSSOVER) {
		delete [] recombination_patterns[0];
		delete [] recombination_patterns;
	}

	recombination_model = FREE_RECOMBINATION;
	if (HG_VERBOSE) cerr<<"...done."<<endl;
	return 0;
}

/**
 * @brief Set a uniform mutation rate for all loci and both directions
 *
 * @param m mutation rate
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_lowd::set_mutation_rates(double m) {
	if (not mem) {
		cerr<<"haploid_lowd::set_mutation_rates(): allocate memory first!"<<endl;
		return HG_MEMERR;	
	}

	for (int fb = 0; fb < 2; fb++)
		for (int locus = 0; locus < number_of_loci; locus++)
			mutation_rates[fb][locus] = m;
	return 0;
}

/**
 * @brief Set two mutation rates (forward / backward) for all loci
 *
 * @param mforward forward mutation rate
 * @param mbackward backward mutation rate
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_lowd::set_mutation_rates(double mforward, double mbackward) {
	if (not mem) {
		cerr<<"haploid_lowd::set_mutation_rates(): allocate memory first!\n";
		return HG_MEMERR;
	}

	for (int locus = 0; locus < number_of_loci; locus++) {
		mutation_rates[0][locus] = mforward;
		mutation_rates[1][locus] = mbackward;
	}
	return 0;
}

/**
 * @brief Set mutation rates (locus specific, both directions the same)
 *
 * @param m array of mutation rates
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_lowd::set_mutation_rates(double* m) {
	if (not mem) {
		cerr<<"haploid_lowd::set_mutation_rates(): allocate memory first!\n";
		return HG_MEMERR;
	}

	for (int locus = 0; locus < number_of_loci; locus++) {
		mutation_rates[0][locus] = m[locus];
		mutation_rates[1][locus] = m[locus];
	}
	return 0;
}

/**
 * @brief Set mutation rates (locus and direction specific)
 *
 * @param m array of mutation rates
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_lowd::set_mutation_rates(double** m) {
	if (not mem) {
		cerr<<"haploid_lowd::set_mutation_rates(): allocate memory first!\n";
		return HG_MEMERR;
	}

	for (int fb = 0; fb < 2; fb++)
		for (int locus = 0; locus < number_of_loci; locus++)
			mutation_rates[fb][locus] = m[fb][locus];
	return 0;
}

/**
 * @brief Set a new recombination model
 *
 * @param rec_model the new recombination model
 *
 * @returns zero if successful, error codes otherwise
 *
 * The allowed models are FREE_RECOMBINATION, SINGLE_CROSSOVER, and CROSSOVERS.
 *
 * Note: this function call involves memory release and allocation.
 */
int haploid_lowd::set_recombination_model(int rec_model) {
	switch(rec_model) {
		case FREE_RECOMBINATION: break;
		case SINGLE_CROSSOVER: break;
		case CROSSOVERS: break;
		default: return HG_BADARG; break;	
	}

	int err = 0;
	if (rec_model != recombination_model) {
		if(recombination_model != FREE_RECOMBINATION)
			err += free_recombination_mem();
		if (rec_model != FREE_RECOMBINATION)
			err += allocate_recombination_mem(rec_model);
		if(err) {
			cerr <<"haploid_lowd::set_recombination_model(): cannot allocate memory for recombination patterns!"<<endl;
			return HG_MEMERR;
		}
	}
	return 0;
}


/**
 * @brief calculate recombination patterns
 *
 * @param rec_rates a vector of recombination rates.
 * @param rec_model an int with the model of recombination to use
 *
 * *Note*: rec_rates must have length L-1 for linear chromosomes, length L for circular ones.
 * *Note*: rec_model must be either CROSSOVERS or, for linear genomes, SINGLE_CROSSOVER
 *
 * @returns zero if successful, error codes otherwise (e.g. out of memory)
 *
 * A routine the calculates the probability of all possible recombination patterns and
 * subpatterns thereof from a vector of recombination rates (rec_rates) passed as argument.
 * It allocated the memory (\f$3^L\f$ or \f$L \: 2^L\f$) and calculates the entire distribution.
 *
 * *Note*: the default value for rec_model is the current recombination model or, if that is FREE_RECOMBINATION,
 * then it's CROSSOVERS.
 */
int haploid_lowd::set_recombination_rates(double *rec_rates, int rec_model) {
	// default value for the recombination model parameter
	if (rec_model == -1) {
		if(recombination_model != FREE_RECOMBINATION)
			rec_model = recombination_model;
		else
			rec_model = CROSSOVERS;
	}

	// one cannot assign recombination rates to free recombination
	if (rec_model == FREE_RECOMBINATION) {
		if(HG_VERBOSE) cerr <<"haploid_lowd::set_recombination_rates(): You cannot set recombination rates to free recombination!"<<endl;
		return HG_BADARG;	
	}

	// recombination model must be in the allowed list of models
	if ((rec_model != CROSSOVERS) and (rec_model != SINGLE_CROSSOVER)) {
		if(HG_VERBOSE) cerr <<"haploid_lowd::set_recombination_rates(): Recombination model not recognized."<<endl;
		return HG_BADARG;	
	}

	// single crossover cannot work for circular genomes, and the crossover probabilities
	// must add up to a number less equal one :-)
	if (rec_model == SINGLE_CROSSOVER) {
		if (circular) {
			if(HG_VERBOSE) cerr <<"haploid_lowd::set_recombination_rates(): Single crossover not available for circular genomes."<<endl;
			return HG_BADARG;	
		}
	
		double sum = 0;
		for (int locus = 0; locus != number_of_loci - 1; locus++)
			sum += rec_rates[locus];
		if(sum > 1) {
			if(HG_VERBOSE) cerr <<"haploid_lowd::set_recombination_rates(): Sum of crossover rates is greater than one!"<<endl;
			return HG_BADARG;	
		}
	}

	// allocate/release memory on changes of recombination model
	double err = set_recombination_model(rec_model);
	if(err) return err;

	//if memory allocation has been successful, calculate the probabilities of recombination based on the model
	if (rec_model == SINGLE_CROSSOVER)
		err = set_recombination_rates_single_crossover(rec_rates);
	else
		err = set_recombination_rates_general(rec_rates);
	return err;
}

/**
 * @brief Set the recombination rates for the single crossover model
 *
 * @param rec_rates array with the recombination rates
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_lowd::set_recombination_rates_single_crossover(double *rec_rates) {
	if(HG_VERBOSE) cerr <<"haploid_lowd::set_recombination_rates_single_crossover()...";

	// recombination_patterns[0][locus] is the probability of the pattern with the crossover immediately AFTER locus
	double sum = 0;
	for (int locus = 0; locus != number_of_loci - 1; locus++) {
		recombination_patterns[0][locus] = 0.5 * rec_rates[locus];
		sum += rec_rates[locus];
	}
	// and recombination_patterns[0][L-1] is the probability of no crossover at all
	recombination_patterns[0][number_of_loci - 1] = 0.5 * (1.0 - sum);

	if(HG_VERBOSE) cerr <<"done."<<endl;
	return 0;
}

/**
 * @brief Set the recombination rates for the general model
 *
 * @param rec_rates array with the recombination rates
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_lowd::set_recombination_rates_general(double *rec_rates) {
	if(HG_VERBOSE) cerr <<"haploid_lowd::set_recombination_rates_general()..."<<endl;

	// The algorithm is divided in two parts:
	// 1. calculate the patterns with all L loci;
	// 2. find patterns with k < L loci via successive marginalizations;
	double * patterns_order_L = recombination_patterns[(1<<number_of_loci) - 1];

	// 1. calculate the probabilities of different crossover realizations
	int strand = 0, newstrand, strandswitches;
	double rr;
	double sum=0;
	for (int i = 0; i < (1<<number_of_loci); i++) {
		patterns_order_L[i]=1.0;

		//parent of the last locus -- equals locus -1 on a circular genome
		strand=(i & (1<<(number_of_loci-1)))>0?1:0;
		strandswitches=0;
		for (int locus = 0; locus < number_of_loci; locus++) {
			newstrand=((i&(1<<locus))>0)?1:0;	//determine the parent of current locus
	
			// For circular genomes all rates are specified.
			// for linear genomes the first rate is set to a large number, e.g. 50
			rr = circular?rec_rates[locus]:(locus?rec_rates[locus - 1]:50);
			if (strand==newstrand)
				patterns_order_L[i] *= (0.5*(1.0 + exp(-2.0*rr)));	//probability of an even number of crossovers
			else {
				patterns_order_L[i] *= (0.5*(1.0 - exp(-2.0*rr)));  //probability of an odd number of crossovers
				strandswitches++;
			}
			strand=newstrand;
		}
		// the constraint of even number of crossovers must be enforced because the genome is
		// always internally represented as circular
		if (strandswitches%2) {patterns_order_L[i] = 0; cerr<<"haploid_lowd::set_recombination_rates_general(): should never happen"<<endl;}
		sum += patterns_order_L[i];
	}
	// normalize the recombination patterns
	if(sum < HG_NOTHING) {
		if(HG_VERBOSE) cerr<<"Recombination rates must be greater than zero!"<<endl;
		return HG_BADARG;
	}

	//normalize the patterns
	for (int i = 0; i < (1<<number_of_loci); i++)
		patterns_order_L[i] /= sum;

	// 2. marginalize the full recombination patterns to the obtain the subsets
	//needed to calculate the recombinant distribution
	int err = marginalize_recombination_patterns();
	if(err) return err;

	if(HG_VERBOSE) cerr <<"done."<<endl;
	return 0;
}

/**
 * @brief set recombination patterns with a list of index value pairs
 *
 * @param vector of index_value_pair_t structures
 *
 * @returns zero if successful, error code otherwise
 */
int haploid_lowd::set_recombination_patterns(vector<index_value_pair_t> iv){
	if(HG_VERBOSE) cerr <<"haploid_lowd::set_recombination_patterns()..."<<endl;

	// allocate/release memory on changes of recombination model
	int err = set_recombination_model(CROSSOVERS);
	if(err) return err;

	// reset the recombination patterns
	double * patterns_order_L = recombination_patterns[(1<<number_of_loci) - 1];
	for (int i = 0; i < (1<<number_of_loci); i++)
		patterns_order_L[i]=0;

	// set the specified recombination patterns
	vector<index_value_pair_t>::iterator pair;
	for (pair=iv.begin();pair!=iv.end(); pair++)
		if((pair->index < (unsigned int)(1<<number_of_loci)) and (pair->val > 0)) {
			//parents are symmetric, hence assign the complementary pattern the same value
			patterns_order_L[pair->index] = pair->val;
			patterns_order_L[~(int)(pair->index)] = pair->val;
		}

	// normalize the recombination patterns
	double sum=0;
	for (int i = 0; i < (1<<number_of_loci); i++)
		sum += patterns_order_L[i];

	for (int i = 0; i < (1<<number_of_loci); i++)
		patterns_order_L[i] /= sum;

	if(sum < HG_NOTHING) {
		if(HG_VERBOSE) cerr<<"Recombination rates must be greater than zero!"<<endl;
		return HG_BADARG;
	}

	// 2. marginalize the full recombination patterns to the obtain the subsets
	//needed to calculate the recombinant distribution
	err+=marginalize_recombination_patterns();
	if(HG_VERBOSE) cerr <<"done."<<endl;
	return err;
}


/**
 * @brief Starting from a fully specified set of probabilities of recombination patterns, compute all sub patterns
 *
 * @param
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_lowd::marginalize_recombination_patterns() {
	if(HG_VERBOSE) cerr <<"haploid_lowd::marginalize_recombination_patterns()..."<<endl;

	int marg_locus, higher_order_subset, higher_order_rec_pattern;
	double *rptemp;

	// marginalize repeatedly until the bottom
	// loop over set of spins of different size, starting with 11111101111 type patterns
	// then 11101110111 type patterns etc. first loop is over different numbers of ones, i.e. spins
	for (int set_size = number_of_loci-1; set_size >= 0; set_size--)
		//loop over all 2^L binary patterns
		for (int subset = 0; subset < (1<<number_of_loci); subset++)
			//if correct number of ones... (its the same in every hypercube...)
			if (fitness.order[subset]==set_size) {
				//determine the first zero, i.e. a locus that can be used to marginalize
				marg_locus=-1;
				for (int locus=0; locus < number_of_loci; locus++)
					if ((subset&(1<<locus))==0)
						{marg_locus=locus; break;}

				//a short hand for the higher order recombination pattern, from which we will marginalize
				higher_order_subset = subset+(1<<marg_locus);
				rptemp = recombination_patterns[higher_order_subset];

				//loop over all pattern of the length set_size and marginalize
				//i.e. 111x01011=111001011+111101011
				for (int rec_pattern = 0; rec_pattern < (1<<set_size); rec_pattern++) {
					higher_order_rec_pattern = (rec_pattern&((1<<marg_locus)-1)) + ((rec_pattern&((1<<set_size)-(1<<marg_locus)))<<1);
					recombination_patterns[subset][rec_pattern] = rptemp[higher_order_rec_pattern] + rptemp[higher_order_rec_pattern+(1<<marg_locus)];
				}
			}

	if(HG_VERBOSE) cerr <<"done."<<endl;
	return 0;
}

/**
 * @brief Get the recombination rate between a locus and the next one.
 *
 * @param locus the locus AFTER which the recombination rate is seeked
 *
 * @returns the recombination rate, -1 for FREE_RECOMBINATION
 *
 * Note: for SINGLE_CROSSOVER, the recombination rate after the last locus is
 *       the rate of no recombination at all.
 */
double haploid_lowd::get_recombination_rate(int locus) {
	if(number_of_loci < 2) {
		if(HG_VERBOSE) cerr<<"There is no recombination with less than 2 loci!"<<endl;
		return HG_BADARG;
	}

	if(recombination_model == FREE_RECOMBINATION)
		return -1;
	else if(recombination_model == SINGLE_CROSSOVER)
		return recombination_patterns[0][locus];
	else if(recombination_model == CROSSOVERS) {
		double *pat = recombination_patterns[(1<<number_of_loci) - 1];
		double r = 0;
		for (int i = 0; i < (1<<number_of_loci); i++)
			if(bool(i&(1<<locus)) xor bool(i&(1<<(locus+1))))
				r += pat[i];
		// Invert the probability of odd crossovers
		r = - 0.5 * log(1 - 2 * r);
		return r;
	}
	else {
		if(HG_VERBOSE) cerr<<"Recombination model not found!"<<endl;
		return HG_BADARG;
	}
}


/**
 * @brief Initialize population in linkage equilibrium.
 *
 * @param freq allele frequencies
 * @param N_in population size
 *
 * @returns zero if successful, error codes otherwise
 *
 * Note: when this function is used to initialize the population, it is likely that the fitness distribution
 * has a very large width. In turn, this can result in an immediate and dramatic drop in diversity within the
 * first few generations. Please check fitness statistics before starting the evolution if this worries you.
 *
 * *Note*: the population size will be set. If not set yet, the carrying
 * capacity will be set to the same number.
 */
int haploid_lowd::set_allele_frequencies(double *freq, unsigned long N_in) {
	population_size=N_in;
	if(carrying_capacity < HG_NOTHING)
		carrying_capacity = population_size;

	// Set the allele frequencies
	population.set_state(HC_FUNC);
	double prob;
	for (int i=0; i<(1<<number_of_loci); i++){
		prob=1.0;
		for (int locus=0; locus < number_of_loci; locus++) {
			if (i&(1<<locus)) prob *= freq[locus];
			else prob *= (1.0-freq[locus]);
		}
		population.func[i] = prob;
	}
	return population.fft_func_to_coeff();
}

/**
 * @brief Initialize the population with genotype counts
 *
 * @param gt vector of index_value_pair_t type with indices and counts
 *
 * @returns zero if successful, error codes otherwise
 *
 * *Note*: the population size will be set as the total sum of counts of all
 * genotypes. If not set yet, the carrying capacity will be set to the same
 * number.
 */
int haploid_lowd::set_genotypes(vector <index_value_pair_t> gt) {
	// Initialize the genotypes
	int err = population.init_list(gt, false);
	if(err) return err;

	// Set the population size as the sum of the clones in input
	population_size=0;
	vector<index_value_pair_t>::iterator pair;
	for (pair = gt.begin(); pair != gt.end(); pair++)
		population_size += pair->val;

	if(carrying_capacity < HG_NOTHING)
		carrying_capacity = population_size;
	return population.normalize();
}

/**
 * @brief Initialize a wildtype population (00...0)
 *
 * @param N_in number of individuals
 *
 * @returns 0 if successful, nonzero otherwise.
 *
 * *Note*: the population size will be set. If not set yet, the carrying
 * capacity will be set to the same number.
 */
int haploid_lowd::set_wildtype(unsigned long N_in) {
	// Initialize the genotype
	index_value_pair_t wildtype(0,N_in);
	vector <index_value_pair_t> temp(1, wildtype);
	int err = population.init_list(temp, false);
	if(err) return err;

	// Set the population size as the sum of the clones in input
	population_size = N_in;
	if(carrying_capacity < HG_NOTHING)
		carrying_capacity = population_size;
	return population.normalize();
}

/**
 * @brief Evolve the population for some generations
 *
 * @param gen number of generations
 *
 * @returns zero if successful, error code in the faulty step otherwise
 *
 * *Note*: the order of selection, mutation, recombination, and resampling could be changed
 * according to needs and beliefs. Note that only recombination calculates the inverse
 * fourier transform of the population. It does so BEFORE the recombination step.
 * To evaluate all allele frequencies and linkage disequilibria, call population.fft_func_to_coeff()
 */
int haploid_lowd::evolve(int gen) {
	if (HG_VERBOSE) cerr<<"haploid_lowd::evolve(int gen)...";
	int err=0, g=0;

	// evolve cycle
	while((err == 0) && (g < gen)) {
		if (HG_VERBOSE) cerr<<"generation "<<generation<<endl;
		if(err==0) err=select();
		if(err==0) err=mutate();
		if(err==0) err=recombine();
		if(err==0) err=resample();
		g++;
		generation++;
		if (generation>HG_LONGTIMEGEN) {generation-=HG_LONGTIMEGEN; long_time_generation+=HG_LONGTIMEGEN;}
	}
	if (HG_VERBOSE) {
		if(err==0) cerr<<"done."<<endl;
		else cerr<<"error "<<err<<"."<<endl;
	}
	return err;
}

/**
 * @brief Evolve the population for some generations, without recombination
 *
 * @param gen number of generations
 *
 * @returns zero if successful, error code in the faulty step otherwise
 *
 * *Note*: the order of selection, mutation, and resampling could be changed
 * according to needs and beliefs.
 */
int haploid_lowd::evolve_norec(int gen) {
	if (HG_VERBOSE) cerr<<"haploid_lowd::evolve_norec(int gen)...";
	int err=0, g=0;

	// evolve cycle
	while((err == 0) && (g < gen)) {
		if (HG_VERBOSE) cerr<<"generation "<<generation<<endl;
		if(err==0) err=select();
		if(err==0) err=mutate();
		if(err==0) err=resample();
		g++;
		generation++;
		if (generation>HG_LONGTIMEGEN) {generation-=HG_LONGTIMEGEN; long_time_generation+=HG_LONGTIMEGEN;}
	}
	if (HG_VERBOSE) {
		if(err==0) cerr<<"done."<<endl;
		else cerr<<"error "<<err<<"."<<endl;
	}
	return err;
}

/**
 * @brief Evolve the population for some generations, without resampling (deterministic)
 *
 * @param gen number of generations
 *
 * @returns zero if successful, error code in the faulty step otherwise
 *
 * *Note*: the order of selection, mutation, and recombination could be changed
 * according to needs and beliefs. Note that only recombination calculates the inverse
 * fourier transform of the population. It does so BEFORE the recombination step.
 * To evaluate all allele frequencies and linkage disequilibria, call population.fft_func_to_coeff()
 */
int haploid_lowd::evolve_deterministic(int gen) {
	if (HG_VERBOSE) cerr<<"haploid_lowd::evolve(int gen)...";
	int err=0, g=0;

	// evolve cycle
	while((err == 0) && (g < gen)) {
		if (HG_VERBOSE) cerr<<"generation "<<generation<<endl;
		if(err==0) err=select();
		if(err==0) err=mutate();
		if(err==0) err=recombine();
		g++;
		generation++;
		if (generation>HG_LONGTIMEGEN) {generation-=HG_LONGTIMEGEN; long_time_generation+=HG_LONGTIMEGEN;}
	}
	if (HG_VERBOSE) {
		if(err==0) cerr<<"done."<<endl;
		else cerr<<"error "<<err<<"."<<endl;
	}
	return err;
}

/**
 * @brief Selection step
 *
 * @returns zero if successful, nonzero otherwise
 *
 * *Note*: Population distribution is reweighted with exp(fitness) and renormalized.
 */
int haploid_lowd::select() {
	population.set_state(HC_FUNC);
	double norm=0;
	for (int i = 0; i < (1<<number_of_loci); i++) {
		population.func[i] *= exp(fitness.func[i]);
		norm += population.func[i];
	}
	population.scale(1.0 / norm);
	return 0;
}

/**
 * @brief Resample the population according to the carrying capacity
 *
 * @returns zero if successful, error codes otherwise
 *
 * *Note*: genotypes with few individuals are sampled using the Poisson distribution, allowing for strict zero;
 * genotypes with many individuals are resampled using a Gaussian distribution, for performance reasons.
 */
int haploid_lowd::resample() {
	population.set_state(HC_FUNC);
	double threshold_HG_CONTINUOUS = double(HG_CONTINUOUS) / carrying_capacity;
	population_size=0;
	//loop over all possible genotypes
	for (int i = 0; i < (1<<number_of_loci); i++) {
		if (population.func[i]<threshold_HG_CONTINUOUS)
			population.func[i] = double(gsl_ran_poisson(rng, carrying_capacity*population.func[i])) / carrying_capacity;
		else
			population.func[i] += double(gsl_ran_gaussian(rng, sqrt(population.func[i]/carrying_capacity)));
		population_size += population.func[i];
	}
	//normalize
	if (population_size<HG_NOTHING)
		return HG_EXTINCT;
	else
		population.scale(1.0 / population_size);
	population_size *= carrying_capacity;
	return 0;
}

/**
 * @brief Mutation step
 *
 * @returns zero if successful, nonzero otherwise
 *
 * Calculate the distribution of mutants and update the population distribution
 */
int haploid_lowd::mutate() {
	mutants.set_state(HC_FUNC);
	population.set_state(HC_FUNC);
	//loop over all possible genotypes
	for (int i = 0; i < (1<<number_of_loci); i++) {
		mutants.func[i] = 0;
		//loop over all possible mutations
		for (int locus=0; locus<number_of_loci; locus++) {
			if (i&(1<<locus))
				mutants.func[i] += mutation_rates[0][locus] * population.func[i-(1<<locus)] - mutation_rates[1][locus] * population.func[i];
			else
				mutants.func[i] += mutation_rates[1][locus] * population.func[i+(1<<locus)] - mutation_rates[0][locus] * population.func[i];
		}
	}
	for (int i = 0; i < (1<<number_of_loci); i++)
		population.func[i] += mutants.func[i];
	return 0;
}

/**
 * @brief Recombination step
 *
 * @returns zero if successful, error codes otherwise
 *
 * Calculate the distribution of recombinants and update the population, a fraction (outcrossing_rate) of the
 * population is replaced by the recombinant distribution, calculated differently depending on whether free
 * recombination or general recombination is used
 */
int haploid_lowd::recombine() {
	int err;
	population.set_state(HC_FUNC);

	if (recombination_model == FREE_RECOMBINATION)
		err=calculate_recombinants_free();
	else if (recombination_model == SINGLE_CROSSOVER)
		err=calculate_recombinants_single();
	else
		err=calculate_recombinants_general();

	for (int i = 0; i < (1<<number_of_loci); i++) {
		population.func[i] += outcrossing_rate * (recombinants.func[i] - population.func[i]);
		// check that genotype frequencies are positive, as
		// negative ones could come from numerical errors in the FFT
		if(population.func[i] < HG_NOTHING)
			population.func[i] = 0;
	}

	return err;
}

/**
 * @brief Calculate the recombinant distribution for the free recombination case
 *
 * @returns zero if successful, nonzero otherwise
 *
 * Almost the same as for the more general case below, but kept separate for
 * performance reasons - this is the most expensive part (3^L).
 */
int haploid_lowd::calculate_recombinants_free() {
	int maternal_alleles, paternal_alleles, count;

	if (HG_VERBOSE) {cerr<<"haploid_lowd::calculate_recombinants_free()...";}

	// prepare hypercubes
	population.fft_func_to_coeff();
	recombinants.set_state(HC_COEFF);

	//normalization of the distribution
	recombinants.coeff[0]=1.0/(1<<number_of_loci);
	if(HG_VERBOSE >= 2) cerr<<0<<"  "<<recombinants.coeff[0]<<endl;

	//loop of all coefficients of the distribution of recombinants
	for (int i = 1; i < (1<<number_of_loci); i++) {
		recombinants.coeff[i]=0;

		//loop over all possible partitions of the loci s1..sk in R^(k)_s1..sk to mother and father
		for (int j = 0; j < (1<<recombinants.order[i]); j++) {
			count=0;
			maternal_alleles=0;
			paternal_alleles=0;

			//build the integers to pull out the maternal and paternal moments
			for (int k = 0; k < number_of_loci; k++)
				if (i&(1<<k)) {
					if (j&(1<<count)) maternal_alleles += (1<<k);
					else paternal_alleles += (1<<k);
					count++;
				}

			//add this particular contribution to the recombinant distribution
			recombinants.coeff[i] += population.coeff[maternal_alleles] * population.coeff[paternal_alleles];
		}

		//normalize: the factor 1<<number_of_loci is due to a peculiarity of the fft algorithm
		recombinants.coeff[i] *= (1<<(number_of_loci - recombinants.order[i]));
	}
	//backtransform to genotype representation
	recombinants.fft_coeff_to_func();

	if (HG_VERBOSE) cerr<<"done."<<endl;
	return 0;
}

/**
 * @brief Calculate the recombinant distribution for the single crossover case
 *
 * @returns zero if successful, nonzero otherwise
 *
 * Calculate the distribution after recombination assumed in random mating with
 * pairs sampled with replacement.
 */
int haploid_lowd::calculate_recombinants_single() {
	if (HG_VERBOSE) cerr<<"haploid_lowd::calculate_recombinants_single()...";

	int maternal_alleles, paternal_alleles, rec_pattern, crossover_point;

	// prepare hypercubes
	population.fft_func_to_coeff();
	recombinants.set_state(HC_COEFF);

	//normalization of the distribution
	recombinants.coeff[0]=1.0/(1<<number_of_loci);
	if(HG_VERBOSE >= 2) cerr<<0<<"  "<<recombinants.coeff[0]<<endl;
	
	double *RP = recombination_patterns[0];

	//loop of all coefficients of the distribution of recombinants
	for (int i=1; i<(1<<number_of_loci); i++) {
		// two things can happen:
		// 1. everything is contributed by the same parent i
		recombinants.coeff[i] = 2*RP[number_of_loci-1] * population.coeff[i] * population.coeff[0];
		if(HG_VERBOSE >= 3) cerr<<i<<"  "<<recombinants.coeff[i]<<"  "<<population.coeff[i]<<endl;
		
		// 2. both parents contribute in variable proportions
		for (crossover_point=0; crossover_point < number_of_loci-1; crossover_point++) {
			rec_pattern=(2<<crossover_point)-1;
			// pick the maternal and paternal Fourier coefficients
			paternal_alleles=(i&(~rec_pattern));
			maternal_alleles=(i&rec_pattern);
			recombinants.coeff[i] += 2*RP[crossover_point] * population.coeff[maternal_alleles] * population.coeff[paternal_alleles];
			if(HG_VERBOSE >= 3) cerr<<i<<"  "<<recombinants.coeff[i]<<"  "<<population.coeff[paternal_alleles]<<endl;
		}		
		
		//normalize: the factor 1<<number_of_loci is due to a peculiarity of the fft algorithm
		recombinants.coeff[i]*=(1<<(number_of_loci));
		if(HG_VERBOSE >= 2) cerr<<i<<"  "<<recombinants.coeff[i]<<endl;
	}
	//backtransform to genotype representation
	recombinants.fft_coeff_to_func();

	if (HG_VERBOSE) cerr<<"done."<<endl;
	return 0;
}

/**
 * @brief Calculate the recombinant distribution for the general case
 *
 * @returns zero if successful, nonzero otherwise
 *
 * Calculate the distribution after recombination assumed in random mating with
 * pairs sampled with replacement.
 */
int haploid_lowd::calculate_recombinants_general() {
	if (HG_VERBOSE) cerr<<"haploid_lowd::calculate_recombinants_general()...";

	int i,j,k, maternal_alleles, paternal_alleles, count;

	// prepare hypercubes
	population.fft_func_to_coeff();
	recombinants.set_state(HC_COEFF);

	//normalization of the distribution
	recombinants.coeff[0]=1.0/(1<<number_of_loci);
	if(HG_VERBOSE >= 2) cerr<<0<<"  "<<recombinants.coeff[0]<<endl;

	//loop of all coefficients of the distribution of recombinants
	for (i=1; i<(1<<number_of_loci); i++) {
		recombinants.coeff[i]=0;

		//loop over all possible partitions of the loci s1..sk in R^(k)_s1..sk to mother and father
		for (j=0; j<(1<<recombinants.order[i]); j++) {
			count=0;
			maternal_alleles=0;
			paternal_alleles=0;

			//build the integers to pull out the maternal and paternal moments
			for (k=0; k<number_of_loci; k++)
				if (i&(1<<k)) {
					if (j&(1<<count)) maternal_alleles+=(1<<k);
					else paternal_alleles+=(1<<k);
					count++;
				}

			//add this particular contribution to the recombinant distribution
			recombinants.coeff[i]+=recombination_patterns[i][j]*population.coeff[maternal_alleles]*population.coeff[paternal_alleles];
			if(HG_VERBOSE >= 3) cerr<<i<<"  "<<recombinants.coeff[i]<<"  "<<population.coeff[paternal_alleles]<<endl;
		}

		//normalize: the factor 1<<number_of_loci is due to a peculiarity of the fft algorithm
		recombinants.coeff[i]*=(1<<(number_of_loci));
		if(HG_VERBOSE >= 2) cerr<<i<<"  "<<recombinants.coeff[i]<<endl;
	}
	//backtransform to genotype representation
	recombinants.fft_coeff_to_func();

	if (HG_VERBOSE) cerr<<"done."<<endl;
	return 0;
}

/**
 * @brief Get the genotype entropy
 *
 * @returns the genotype entropy of the population
 *
 * The genotype entropy is defined as follows:
 * \f[ S := - \sum_{g} \nu_g \cdot \log \nu_g, \f]
 * where \f$g\f$ runs over all possible genomes, and \f$\nu_g\f$ is the frequency of that genome
 * in the population.
 */
double haploid_lowd::genotype_entropy() {
	// make sure the population is in the right state
	if (population.get_state()==HC_COEFF) population.fft_coeff_to_func();

	double entropy = 0;
	for (int i = 0; i < (1<<number_of_loci); i++)
		entropy -= population.func[i] * log(population.func[i]);
	return entropy;
}

/**
 * @brief Get the allele entropy
 *
 * @returns the allele entropy of the population
 *
 * The allele entropy is defined as follows:
 * \f[ S := - \sum_{i=1}^L \nu_i\log \nu_i + (1 - \nu_i)\log (1 - \nu_i), \f]
 * where \f$\nu_i\f$ is the frequency of the i-th allele.
 */
double haploid_lowd::allele_entropy() {
	// make sure the population is in the right state
	if (population.get_state()==HC_FUNC) population.fft_func_to_coeff();

	double entropy = 0;
	for (int locus = 0; locus < number_of_loci; locus++) {
		entropy -= 0.5 * (1.0 + population.coeff[(1<<locus)]) * log(0.5*(1.0 + population.coeff[(1<<locus)]));
		entropy -= 0.5 * (1.0 - population.coeff[(1<<locus)]) * log(0.5*(1.0 - population.coeff[(1<<locus)]));
	}
	return entropy;
}

/**
 * @brief Get fitness mean and variance in the population
 *
 * @returns stat_t with the requested statistics
 */
stat_t haploid_lowd::get_fitness_statistics() {
	if (population.get_state()==HC_COEFF) population.fft_coeff_to_func();

	double mf = 0, sq = 0, temp;
	for (int i = 0; i < (1<<number_of_loci); i++) {
		temp = population.get_func(i) * fitness.get_func(i);
		mf += temp;
		sq += temp*temp;
	}
	return stat_t(mf, sq - mf * mf);
}
