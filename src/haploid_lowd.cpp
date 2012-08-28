/*
 * haploid_lowd.cpp
 *
 *  Created on: Jan 27, 2010
 *  Author: Richard Neher
 *  Modified by: Fabio Zanini
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
haploid_lowd::haploid_lowd(int L_in, int rng_seed) : number_of_loci(L_in), population_size(0), mem(false), recombination_mem(FREE_RECOMBINATION), recombination_model(FREE_RECOMBINATION), outcrossing_rate(1.0), circular(false), generation(0), long_time_generation(0) {
	if (L_in <1) {
		cerr <<"haploid_lowd::haploid_lowd(): Bad Arguments! L must be larger or equal one."<<endl;
		throw HG_BADARG;
        }

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

	int err=0;
	rng=gsl_rng_alloc(RNG);
	gsl_rng_set(rng, seed);
	cerr <<"haploid_lowd() random number seed: "<<seed<<endl;
	err+=fitness.set_up(number_of_loci, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
	err+=population.set_up(number_of_loci, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
	err+=mutants.set_up(number_of_loci, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
	err+=recombinants.set_up(number_of_loci, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
	mutation_rates=new double*[2]; //allocate backward and forward mutation rates arrays
	mutation_rates[0]=new double [number_of_loci]; //allocate forward mutation rates
	mutation_rates[1]=new double [number_of_loci]; //allocate backward mutation rates
	for (int fb=0; fb<2; fb++){	//set mutation rates to zero
		for (int locus=0; locus<number_of_loci; locus++){
			mutation_rates[fb][locus]=0;
		}
	}
	if (err==0) {
		mem=true;
		return 0;
	}
	else return HG_MEMERR;
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
 * *Note*: this function also sets the recombination_mem flag.
 *
 * *Note*: in the SINGLE_CROSSOVER model, the recombination patterns are sorted as follows:
 *
 *    patterns[0] = 0..000, patterns[L]   = 1..111
 *    patterns[1] = 0..001, patterns[L+1] = 1..110
 *    patterns[2] = 0..011, patterns[L+2] = 1..100
 *    [...]
 */
int haploid_lowd::allocate_recombination_mem(int rec_model) {
		int err = 0;
		int i, spin;
		int temp;
		int *nspins;	//temporary variables the track the number of ones in the binary representation of i
		nspins=new int [1<<number_of_loci];
		if (nspins==NULL) {
			cerr<<"haploid_lowd::set_recombination_rates(): Can not allocate memory!"<<endl;
			return HG_MEMERR;
		}
		spin=-1;
		nspins[0]=0;
		//allocate space for all possible subsets of loci
		recombination_patterns=new double* [1<<number_of_loci];
		recombination_patterns[0]=new double [1];
		//loop over all possible locus subsets and allocate space for all
		//possible ways to assign the subset to father and mother (2^nspins)
		for (i=1; i<(1<<number_of_loci); i++){

			// coefficients are sorted in a specific order
			// the order of coefficient k is 1+(the order of coefficient[k-2^spin])
			if (i==(1<<(spin+1))) spin++;
			temp=1+nspins[i-(1<<spin)];
			nspins[i]=temp;

			//all possible ways to assign the subset to father and mother
			// for CROSSOVERS: 2^temp
			//     e.g. 000, 001, 010, 011, 100, 101, 110, 111
			//
			// for SINGLE_CROSSOVER: 2 * temp
			//     e.g. 000, 001, 011, 111, 110, 100
			if (rec_model == CROSSOVERS)
				recombination_patterns[i]=new double [(1<<temp)];
			else
				recombination_patterns[i]=new double [2 * temp];

			if (recombination_patterns[i]==NULL) err+=1;
		}
		delete [] nspins;
		recombination_mem = CROSSOVERS;

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
	if (recombination_mem != FREE_RECOMBINATION) {
		for (int i=0; i<(1<<number_of_loci); i++){
			delete [] recombination_patterns[i];
		}
		delete [] recombination_patterns;
	}
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
	if (mem){
		for (int fb=0; fb<2; fb++){
			for (int locus=0; locus<number_of_loci; locus++){
				mutation_rates[fb][locus]=m;
			}
		}
		return 0;
	} else {
		cerr<<"haploid_lowd::set_mutation_rates(): allocate memory first!\n";
		return HG_MEMERR;
	}
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
	if (mem){
		for (int locus=0; locus<number_of_loci; locus++){
			mutation_rates[0][locus]=mforward;
			mutation_rates[1][locus]=mbackward;
		}
		return 0;
	} else {
		cerr<<"haploid_lowd::set_mutation_rates(): allocate memory first!\n";
		return HG_MEMERR;
	}
}

/**
 * @brief Set mutation rates (locus specific, both directions the same)
 *
 * @param m array of mutation rates
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_lowd::set_mutation_rates(double* m) {
	if (mem){
		for (int locus=0; locus<number_of_loci; locus++){
			mutation_rates[0][locus]=m[locus];
			mutation_rates[1][locus]=m[locus];
		}
		return 0;
	}else{
		cerr<<"haploid_lowd::set_mutation_rates(): allocate memory first!\n";
		return HG_MEMERR;
	}
}

/**
 * @brief Set mutation rates (locus and direction specific)
 *
 * @param m array of mutation rates
 *
 * @returns zero if successful, error codes otherwise
 */
int haploid_lowd::set_mutation_rates(double** m) {
	if (mem){
		for (int fb=0; fb<2; fb++){
			for (int locus=0; locus<number_of_loci; locus++){
				mutation_rates[fb][locus]=m[fb][locus];
			}
		}
		return 0;
	}else{
		cerr<<"haploid_lowd::set_mutation_rates(): allocate memory first!\n";
		return HG_MEMERR;
	}
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
 * It allocated the memory (\f$3^L\f$) and calculates the entire distribution.
 */
int haploid_lowd::set_recombination_rates(double *rec_rates, int rec_model) {
	double err=0;
	double sum;

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

	if (circular and (rec_model == SINGLE_CROSSOVER)) {
		if(HG_VERBOSE) cerr <<"haploid_lowd::set_recombination_rates(): Single crossover not available for circular genomes."<<endl;
		return HG_BADARG;	
	}

	// check that the sum of probablity densities is not larger than one
	if (rec_model == SINGLE_CROSSOVER) {
		sum = 0;
		for (double * rtmp = rec_rates; rtmp != rec_rates + number_of_loci - 1; rtmp++)
			sum += *rtmp;
		if(sum > 1) {
			if(HG_VERBOSE) cerr <<"haploid_lowd::set_recombination_rates(): Rate of NO crossover less than zero!"<<endl;
			return HG_BADARG;	
		}
	}

	// allocate/release memory on changes of recombination model
	if (rec_model != recombination_mem) {
		if (recombination_mem != FREE_RECOMBINATION)
			err += free_recombination_mem();
		err += allocate_recombination_mem(rec_model);
		// If memory had problems, exit
		if(err) {
			cerr <<"haploid_lowd::set_recombination_rates(): cannot allocate memory for recombination patterns!"<<endl;
			return HG_MEMERR;
		}
	}

	//if memory allocation has been successful, calculate the probabilities of recombination
	// The algorithm is divided in two parts:
	// 1. calculate the patterns with all L loci;
	// 2. find patterns with k < L loci via successive marginalizations;
	int locus;
	double * patterns_order_L = recombination_patterns[(1<<number_of_loci) - 1];
	int subset, marg_locus, higher_order_subset, higher_order_rec_pattern;
	double *rptemp;

	if (rec_model == CROSSOVERS) {
		// 1. calculate the probabilities of different crossover realizations
		int strand=0, newstrand, strandswitches;
		double rr;
		sum=0;
		for (int i=0; i < (1<<number_of_loci); i++) {
			patterns_order_L[i]=1.0;
			strand=(i & (1<<(number_of_loci-1)))>0?1:0;
			strandswitches=0;
			for (locus=0; locus < number_of_loci; locus++) {
				newstrand=((i&(1<<locus))>0)?1:0;
	
				// Circular genomes have all rates, linear ones lack the first (which must be a large number, e.g. 50)
				rr = circular?rec_rates[locus]:(locus?rec_rates[locus - 1]:50);
				if (strand==newstrand)
					patterns_order_L[i] *= (0.5*(1.0 + exp(-2.0*rr)));
				else {
					patterns_order_L[i] *= (0.5*(1.0 - exp(-2.0*rr)));
					strandswitches++;
				}
				strand=newstrand;
			}
			// the constraint of even number of crossovers must be enforced because the genome is
			// always internally represented as circular
			if (strandswitches%2) patterns_order_L[i] = 0;
			sum += patterns_order_L[i];
		}
		for (int i=0; i < (1<<number_of_loci); i++)
			patterns_order_L[i]/=sum;
	
		// 2. marginalize repeatedly until the bottom
		//loop over set of spins of different size, starting with 11111101111 type patterns
		//then 11101110111 type patterns etc. first loop is over different numbers of ones, i.e. spins
		for (int set_size=number_of_loci-1; set_size>=0; set_size--) {
			//loop over all 2^L binary patterns
			for (subset=0; subset < (1<<number_of_loci); subset++) {
				//if correct number of ones... (its the same in every hypercube...)
				if (fitness.order[subset]==set_size) {
					//determine the first zero, i.e. a locus that can be used to marginalize
					marg_locus=-1;
					for (locus=0; locus < number_of_loci; locus++) {
						if ((subset&(1<<locus))==0)
							{marg_locus=locus; break;}
					}
					//a short hand for the higher order recombination pattern, from which we will marginalize
					higher_order_subset=subset+(1<<marg_locus);
					rptemp=recombination_patterns[higher_order_subset];
					//loop over all pattern of the length set_size and marginalize
					//i.e. 111x01011=111001011+111101011
					for (int rec_pattern=0; rec_pattern<(1<<set_size); rec_pattern++) {
						higher_order_rec_pattern=(rec_pattern&((1<<marg_locus)-1))+((rec_pattern&((1<<set_size)-(1<<marg_locus)))<<1);
						recombination_patterns[subset][rec_pattern]=rptemp[higher_order_rec_pattern]+rptemp[higher_order_rec_pattern+(1<<marg_locus)];
					}
				}
			}
		}
		recombination_model = CROSSOVERS;
	}

	// single crossover for linear genomes
	else {
		vector <int> ii;

		// 1. calculate the probabilities of different crossover realizations
		// pat[locus] is the proability of the pattern with the crossover
		// immediately AFTER locus, and pat[L-1] is the probability of
		// no crossover at all
		patterns_order_L[number_of_loci - 1] = patterns_order_L[2 * number_of_loci - 1] = 0.5 * (1.0 - sum);
		for (locus=0; locus < number_of_loci - 1; locus++)
			patterns_order_L[locus] = patterns_order_L[locus + number_of_loci] = 0.5 * rec_rates[locus];

		// 2. marginalize repeatedly until the bottom
		for (int set_size=number_of_loci-1; set_size > 0; set_size--) {
			// recombination_patterns of set_size have length 2 * set_size
			// we have to do three things:
			// I. figure out what set of sites correspond to i
			// II. run over all possible crossover points in the reduced set
			// III. marginalize over all possible invisible crossover points (1 or 2)
			//loop over all 2^L binary patterns
			for (subset=0; subset < (1<<number_of_loci); subset++) {
				//if correct number of ones... (its the same in every hypercube...)
				if (fitness.order[subset] == set_size) {
					//determine the first zero, i.e. a locus that can be used to marginalize
					marg_locus=-1;
					for (locus=0; locus < number_of_loci; locus++) {
						if ((subset&(1<<locus))==0)
							{marg_locus=locus; break;}
					}
					//a short hand for the higher order recombination pattern, from which we will marginalize
					higher_order_subset=subset+(1<<marg_locus);
					rptemp=recombination_patterns[higher_order_subset];

					// I. create the multiindex i = indices(subset != 0)
					ii.clear();
					for (locus=0; locus < number_of_loci; locus++) {
						if ((subset&(1<<locus)))
							ii.push_back(locus);
					}

					// II. loop over all pattern of the length set_size and marginalize
					// This can only happen at the boundary between 0s and 1s
					// i.e. 111x00 = 111000 + 111100
					// but  110x00 = 110000 by virtue of the single crossover!
					// for locus < set_size - 1, the restricted pattern is heterogeneous,
					// e.g. x00xx1x
					for (locus = 0; locus < set_size - 1; locus++) {

						// III. only zero OK
						if (marg_locus < ii[locus])
							recombination_patterns[subset][locus] = rptemp[locus + 1];

						// III. only one OK
						else if (marg_locus > ii[locus + 1])
							recombination_patterns[subset][locus] = rptemp[locus];
						
						// III. both OK
						else
							recombination_patterns[subset][locus] = rptemp[locus + 1] + rptemp[locus];
						
						// there is symmetry between the probabilities 0x0xx1x and 1x1xx0x
						recombination_patterns[subset][locus + set_size] = recombination_patterns[subset][locus];
					}
					// II. for locus == set_size - 1, the restricted pattern is homogeneous,
					// e.g. x00xx0x
					// III. both OK
					if ((marg_locus < ii[0]) or (marg_locus > ii[set_size - 1]))
						recombination_patterns[subset][set_size - 1] = rptemp[set_size] + rptemp[2 * set_size + 1];
					// III. only zero OK
					else
						recombination_patterns[subset][set_size - 1] = rptemp[set_size];
					
					// there is symmetry between the probabilities 0x0xx1x and 1x1xx0x
					recombination_patterns[subset][2 * set_size - 1] = recombination_patterns[subset][set_size - 1];
				}			
			}
		}
		// the very bottom has len(ii) = 0, hence deserves a special treatment
		recombination_patterns[0][0] = recombination_patterns[subset][0] + recombination_patterns[subset][1];
		recombination_model = SINGLE_CROSSOVER;
	}
	return 0;
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
        // Set the population size
	population_size=N_in;
        if(carrying_capacity < HG_NOTHING)
                carrying_capacity = population_size;

        // Set the allele frequencies
	int locus, i;
	double prob;
	population.set_state(HC_FUNC);
	for (i=0; i<(1<<number_of_loci); i++){
		prob=1.0;
		for (locus=0; locus<number_of_loci; locus++){
			if (i&(1<<locus)) prob*=freq[locus];
			else prob*=(1.0-freq[locus]);
		}
		population.func[i]=prob;
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
	for(size_t i = 0; i < gt.size(); i++)
		population_size += gt[i].val;
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

        // Initialize the genotypes
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
	for (int i=0; i<(1<<number_of_loci); i++){
		population.func[i]*=exp(fitness.func[i]);
		norm+=population.func[i];
	}
	population.scale(1.0/norm);
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
	double threshold_HG_CONTINUOUS=double(HG_CONTINUOUS)/carrying_capacity;
	population_size=0;
	for (int i=0; i<(1<<number_of_loci); i++){
		if (population.func[i]<threshold_HG_CONTINUOUS)
		{
			population.func[i]=double(gsl_ran_poisson(rng, carrying_capacity*population.func[i]))/carrying_capacity;
		}
		else
		{
			population.func[i]+=double(gsl_ran_gaussian(rng, sqrt(population.func[i]/carrying_capacity)));
		}
		population_size += population.func[i];
	}
	if (population_size<HG_NOTHING){
		return HG_EXTINCT;
	}
	else population.scale(1.0/population_size);
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
	int locus;
	mutants.set_state(HC_FUNC);
	population.set_state(HC_FUNC);
	for (int i=0; i<(1<<number_of_loci); i++) {
		mutants.func[i]=0;
		for (locus=0; locus<number_of_loci; locus++) {
			if (i&(1<<locus))
				mutants.func[i]+=mutation_rates[0][locus]*population.func[i-(1<<locus)]-mutation_rates[1][locus]*population.func[i];
			else
				mutants.func[i]+=mutation_rates[1][locus]*population.func[i+(1<<locus)]-mutation_rates[0][locus]*population.func[i];
		}
	}
	for (int i=0; i<(1<<number_of_loci); i++){
		population.func[i]+=mutants.func[i];
	}
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
	if (recombination_model == FREE_RECOMBINATION){
		err=calculate_recombinants_free();
		for (int i=0; i<(1<<number_of_loci); i++){
			population.func[i]+=outcrossing_rate*(recombinants.func[i]-population.func[i]);
		}
	}else{
		err=calculate_recombinants_general();
		for (int i=0; i<(1<<number_of_loci); i++){
			population.func[i]+=outcrossing_rate*(recombinants.func[i]-population.func[i]);
		}
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
	int i,j,k, maternal_alleles, paternal_alleles, count;

	// prepare hypercubes
	population.fft_func_to_coeff();
	recombinants.set_state(HC_COEFF);

	//normalization of the distribution
	recombinants.coeff[0]=1.0/(1<<number_of_loci);

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
			recombinants.coeff[i]+=population.coeff[maternal_alleles]*population.coeff[paternal_alleles];
		}

		//normalize: the factor 1<<number_of_loci is due to a peculiarity of the fft algorithm
		recombinants.coeff[i]*=1.0*(1<<(number_of_loci-recombinants.order[i]));
	}

	//backtransform to genotype representation
	recombinants.fft_coeff_to_func();
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
	int i,j,k, maternal_alleles, paternal_alleles, count;

	// prepare hypercubes
	population.fft_func_to_coeff();
	recombinants.set_state(HC_COEFF);

	//normalization of the distribution
	recombinants.coeff[0]=1.0/(1<<number_of_loci);
	if(HG_VERBOSE) cerr<<0<<"  "<<recombinants.coeff[0]<<endl;

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
			if(HG_VERBOSE >= 2) cerr<<i<<"  "<<recombinants.coeff[i]<<"  "<<population.coeff[paternal_alleles]<<endl;
		}

		//normalize: the factor 1<<number_of_loci is due to a peculiarity of the fft algorithm
		recombinants.coeff[i]*=(1<<(number_of_loci));
		if(HG_VERBOSE) cerr<<i<<"  "<<recombinants.coeff[i]<<endl;
	}

	//backtransform to genotype representation
	recombinants.fft_coeff_to_func();
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
double haploid_lowd::genotype_entropy(){
	double S=0;
	if (population.get_state()==HC_COEFF) population.fft_coeff_to_func();
	for (int i=0; i<(1<<number_of_loci); i++){
		S-=population.func[i]*log(population.func[i]);
	}
	return S;
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
double haploid_lowd::allele_entropy(){
	double SA=0;
	if (population.get_state()==HC_FUNC) population.fft_func_to_coeff();
	for (int locus=0; locus<number_of_loci; locus++){
		SA-=0.5*(1.0+population.coeff[(1<<locus)])*log(0.5*(1.0+population.coeff[(1<<locus)]));
		SA-=0.5*(1.0-population.coeff[(1<<locus)])*log(0.5*(1.0-population.coeff[(1<<locus)]));
	}
	return SA;
}

/**
 * @brief Get fitness mean and variance in the population
 *
 * @returns stat_t with the requested statistics
 */
stat_t haploid_lowd::get_fitness_statistics(){
	double mf=0, sq=0, temp;
	if (population.get_state()==HC_COEFF) population.fft_coeff_to_func();
	for (int locus=0; locus<1<<number_of_loci; locus++){
		temp=population.get_func(locus)*fitness.get_func(locus);
		mf+=temp;
		sq+=temp*temp;
	}
	return stat_t(mf, sq-mf);
}


/**
 * @brief Test the recombination routine using Fourier transforms
 *
 * @returns zero if both routines agree, -1 otherwise
 *
 * Debugging routine: calculates the distribution of recombinants explicitly and
 * compares the result to the recombinant distribution obtained via fourier transform
 */
int haploid_lowd_test::test_recombinant_distribution(){
	double *test_rec;
	double dev=0;
	//allocate memory for the recombinant distribution calculated step-by-step
	test_rec=new double [(1<<number_of_loci)];
	int mother, father;
	//now calculate the recombinant distribution from pairs of parents.
	int gt1, gt2, rec_pattern;
	if (recombination_model == FREE_RECOMBINATION){
		//calculate recombinants the efficient way
		calculate_recombinants_free();
		for (gt1=0; gt1<(1<<number_of_loci); gt1++){	//target genotype
			test_rec[gt1]=0.0;							//initialize
			//loop over all recombination patterns (equal probability)
			for (rec_pattern=0; rec_pattern<(1<<number_of_loci); rec_pattern++){
				//loop over the parts of the maternal and paternal genomes not inherited
				for (gt2=0; gt2<(1<<number_of_loci); gt2++){
					//construct maternal and paternal genotypes
					mother=(gt1&(rec_pattern))+(gt2&(~rec_pattern));
					father=(gt1&(~rec_pattern))+(gt2&(rec_pattern));
					//increment the rec distribution
					test_rec[gt1]+=population.func[mother]*population.func[father];
				}
			}
			//normalize
			test_rec[gt1]*=1.0/(1<<number_of_loci);
			cout <<gt1<<"  "<<test_rec[gt1]<<"  "<<recombinants.func[gt1]<<endl;
			//sum up all deviations
			dev+=(test_rec[gt1]-recombinants.func[gt1])*(test_rec[gt1]-recombinants.func[gt1]);
		}
	}else{ 	//same as above, only the individual contribution
		//calculate recombinants the efficient way
		calculate_recombinants_general();
		for (gt1=0; gt1<(1<<number_of_loci); gt1++){
			test_rec[gt1]=0.0;
			for (rec_pattern=0; rec_pattern<(1<<number_of_loci); rec_pattern++){
				for (gt2=0; gt2<(1<<number_of_loci); gt2++){
					mother=(gt1&(rec_pattern))+(gt2&(~rec_pattern));
					father=(gt1&(~rec_pattern))+(gt2&(rec_pattern));
					//contribution is weighted by the probability of this particular recombination pattern
					//this got calculated and stored in recombination_patterns[(1<<number_of_loci)-1]
					test_rec[gt1]+=recombination_patterns[(1<<number_of_loci)-1][rec_pattern]*population.func[mother]*population.func[father];
				}
			}
			cout <<gt1<<"  "<<test_rec[gt1]<<"  "<<recombinants.func[gt1]<<endl;
			dev+=(test_rec[gt1]-recombinants.func[gt1])*(test_rec[gt1]-recombinants.func[gt1]);
		}
	}
	delete [] test_rec;
	if (dev>1e-9){
		cout <<"Deviation between explicit and fourier transform version! "<<dev<<endl;
		return -1;
	}else{
		cout <<"Explicit and fourier transform version agree to "<<dev<<endl;
		return 0;
	}
	return 0;
}

/**
 * @brief Test the recombination routine extensively
 *
 * @param rec_rates recombination rates used for testing
 *
 * @returns zero (but look at the stdout)
 *
 * Debugging routine: produces random genotypes configurations and test whether they recombine correctly.
 */
int haploid_lowd_test::test_recombination(double *rec_rates){

	//calculate the genetic map, i.e. cumulative recombination rates
	double* cumulative_rates=new double [number_of_loci+1];
	cumulative_rates[0]=0.0;
	for (int locus =1; locus<number_of_loci+1; locus++) cumulative_rates[locus]=cumulative_rates[locus-1]+rec_rates[locus-1];

	//initialize the internal recombination rates
	set_recombination_rates(rec_rates);

	//initialize the population randomly and test the recombination procedure
	population.set_state(HC_FUNC);
	for (int r=0; r<1; r++){
		for (int i=0; i<(1<<number_of_loci); i++){
			population.func[i]=gsl_rng_uniform(rng);
		}
		population.normalize();
		test_recombinant_distribution();
	}
	population.set_state(HC_FUNC);
	for (int i=0; i<(1<<number_of_loci); i++){
		population.func[i]=gsl_rng_uniform(rng);
	}
	population.normalize();

	//study the decay of cumulants from the randomly initialized initialized population
	//output header
	cout <<"\n\nRatio of the cumulants and the expected decay curve, should be constant. Last column shows dynamic range\n";
	cout <<"Generation  ";
	for (int l1=0; l1<number_of_loci; l1++){
		for(int l2=0; l2<l1; l2++){
			cout <<setw(13)<<l1<<" "<<l2;
		}
	}
	cout <<setw(15)<<"exp(-rmax*t)";
	cout<<'\n';
	//for a thousand time steps, recombine and watch the cumulants decay
	for (int g=0; g<1000; g++){
		if (g%100==0){ //output every hundred generations
			cout <<setw(10)<<g;
			for (int l1=0; l1<number_of_loci; l1++){
				for(int l2=0; l2<l1; l2++){
					cout <<setw(15)<<get_LD(l1,l2)*exp(g*0.5*(1.0-exp(-2.0*(cumulative_rates[l1+1]-cumulative_rates[l2+1]))));
				}
			}
			cout <<setw(15)<<exp(-g*0.5*(1.0-exp(-2.0*(cumulative_rates[number_of_loci]-cumulative_rates[1]))));
			cout<<'\n';
		}
		recombine();
	}
	delete cumulative_rates;
	return 0;
}



/**
 * @brief Test the mutation-drift equilibrium with diffusion theory
 *
 * @param mu mutation rates
 *
 * @returns zero (but look at the stdout)
 */
int haploid_lowd_test::mutation_drift_equilibrium(double **mu){
	set_mutation_rates(mu);

	//init population and recombination rates
	double *af=new double[number_of_loci];;
	double *recrates=new double[number_of_loci];;
	for (int i=0; i<number_of_loci; i++){
		af[i]=0;
		recrates[i]=10;
	}
	set_allele_frequencies(af, 1000);
	//allocate histograms to store allele frequency distributions
	gsl_histogram **mutfreq=new gsl_histogram* [number_of_loci];
	for (int locus=0; locus<number_of_loci; locus++){
		mutfreq[locus]=gsl_histogram_alloc(100);
		gsl_histogram_set_ranges_uniform(mutfreq[locus], -1,1);
	}

	//equilibrate for 2N generations
	for (int gen=0; gen<2*carrying_capacity; gen++){
		mutate();
		resample();
	}
	//take 100000 samples every 1000 generations (assumes population is of order 1000)
	for (int r=0; r<100000; r++){
		for (int gen=0; gen<1000; gen++){
			mutate();
			resample();
		}

		for (int locus=0; locus<number_of_loci; locus++){
			gsl_histogram_increment(mutfreq[locus], get_chi(locus));
		}
	}

	//output: normalized histograms as well as theoretical expectation from diffusion theory
	//calculate norm of distributions first, output below.
	double upper, lower;
	double* histogramnorm=new double [number_of_loci];
	double* theorynorm=new double [number_of_loci];
	for (int locus=0; locus<number_of_loci; locus++){
		histogramnorm[locus]=0;
		theorynorm[locus]=0;
		for (int i=0; i<100; i++){
			gsl_histogram_get_range(mutfreq[locus], i, &lower, &upper);
			histogramnorm[locus]+=gsl_histogram_get(mutfreq[locus], i);
			theorynorm[locus]+=pow(0.5*(1+0.5*(upper+lower)), 2*carrying_capacity*mu[0][locus]-1)*pow(0.5*(1-0.5*(upper+lower)), 2*carrying_capacity*mu[1][locus]-1);
		}
	}
	for (int i=0; i<100; i++){
		gsl_histogram_get_range(mutfreq[0], i, &lower, &upper);
		cout <<setw(15)<<0.5*(upper+lower);
		for (int locus=0; locus<number_of_loci; locus++){
			cout <<setw(15)<<gsl_histogram_get(mutfreq[locus], i)/histogramnorm[locus]
					<<setw(15)<<pow(0.5*(1+0.5*(upper+lower)), 2*carrying_capacity*mu[0][locus]-1)*pow(0.5*(1-0.5*(upper+lower)), 2*carrying_capacity*mu[1][locus]-1)/theorynorm[locus];
		}
		cout <<endl;
	}
	return 0;
}
