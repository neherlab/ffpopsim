/*
 * haploid_gt_dis.cpp
 *
 *  Created on: Jan 27, 2010
 *  Author: Richard Neher & Boris Shraiman
 */
#include "popgen_lowd.h"

haploid_gt_dis::haploid_gt_dis()
{
	number_of_loci=0;
	free_recombination=true;
	outcrossing_rate=0.0;
	generation=0;
	long_time_generation=0.0;
	circular=false;
	mem=false;
}

haploid_gt_dis::~haploid_gt_dis()
{
	free_mem();
}

haploid_gt_dis::haploid_gt_dis(double N_in, int L_in, int rngseed)
{
	mem=false;
	free_recombination=true;
	outcrossing_rate=0.0;
	generation=0;
	long_time_generation=0.0;
	circular=false;
	set_up(N_in, L_in, rngseed);
}

int haploid_gt_dis::set_up(double N_in, int L_in, int rngseed)
{
	population_size=N_in;
	number_of_loci=L_in;
	if (rngseed==0) seed=time(NULL);
	else seed=rngseed;
	return allocate_mem();
}

/*
 * allocate memory and set up the different hypercubes needed to store the fitness, population
 * recombinants, and mutants
 */
int haploid_gt_dis::allocate_mem()
{
	int err=0;
	rng=gsl_rng_alloc(RNG);
	gsl_rng_set(rng, seed);
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


/*
 * Free memory allocated in allocate_mem and set_recombination_rates
 */
int haploid_gt_dis::free_mem()
{
	if (mem) {
		fitness.~hypercube();
		population.~hypercube();
		recombinants.~hypercube();
		mutants.~hypercube();
		gsl_rng_free(rng);
		if (!free_recombination){
			for (int i=0; i<(1<<number_of_loci); i++){
				delete [] recombination_patters[i];
			}
			delete [] recombination_patters;
		}
		delete [] mutation_rates[0];
		delete [] mutation_rates[1];
		delete [] mutation_rates;
	}
	mem=false;
	return 0;
}

/*
 * Initialize the population in linkage equilibrium
 */
int haploid_gt_dis::init_frequencies(double *freq){
	double prob;
	int locus, i;
	population.set_state(HC_FUNC);
	for (i=0; i<(1<<number_of_loci); i++){
		prob=1.0;
		for (locus=0; locus<number_of_loci; locus++){
			if (i&(1<<locus)) prob*=freq[locus];
			else prob*=(1.0-freq[locus]);
		}
		population.func[i]=prob;
	}
	generation=0;
	long_time_generation=0.0;
	return population.fft_func_to_coeff();
}

/*
 * Initialize the population in linkage equilibrium
 */
int haploid_gt_dis::init_genotypes(vector <index_value_pair_t> gt){
	population.init_list(gt, false);
	generation=0;
	long_time_generation=0.0;
	return population.normalize();
}

/**
 * @brief evolve the population for gen generations
 *
 * Evolve the population for one generation, for finite and infinite populations
 * The order of selection, mutation, recombination, and resampling could be changed
 * according to needs and beliefs. Note that only recombination calculates the inverse
 * fourier transform of the population. It does so BEFORE the recombination step.
 * To evaluate all allele frequencies and linkage disequilibria, call population.fft_func_to_coeff()
 *
 * @param gen number of generations
 *
 * @returns sum of error codes for the four steps (selection, mutation, recombination, resampling)
 */
int haploid_gt_dis::evolve(int gen){
	int err=0;
	for (int g=0; g<gen; g++){
		err+=select();
		err+=mutate();
		err+=recombine();
		err+=resample();
	}
	generation+=gen;
	if (generation>HG_LONGTIMEGEN) {generation-=HG_LONGTIMEGEN; long_time_generation+=HG_LONGTIMEGEN;}
	return err;
}

int haploid_gt_dis::evolve_norec(int gen){
	int err=0;
	for (int g=0; g<gen; g++){
		err+=select();
		err+=mutate();
		err+=resample();
	}
	generation+=gen;
	if (generation>HG_LONGTIMEGEN) {generation-=HG_LONGTIMEGEN; long_time_generation+=HG_LONGTIMEGEN;}
	return err;
}

int haploid_gt_dis::evolve_deterministic(int gen){
	int err=0;
	for (int g=0; g<gen; g++){
		err+=select();
		err+=mutate();
		err+=recombine();
	}
	generation+=gen;
	if (generation>HG_LONGTIMEGEN) {generation-=HG_LONGTIMEGEN; long_time_generation+=HG_LONGTIMEGEN;}
	return err;
}

/*
 * Selection step: Population distribution is reweighted with exp(fitness) and renormalized
 */
int haploid_gt_dis::select()
{
	population.set_state(HC_FUNC);
	double norm=0;
	for (int i=0; i<(1<<number_of_loci); i++){
		population.func[i]*=exp(fitness.func[i]);
		norm+=population.func[i];
	}
	population.scale(1.0/norm);
	return 0;
}

/*
 * resample the population distribution to produce a population of approx N discrete individuals
 * genotypes with few individuals are sampled using the Poisson distribution, allowing for strict zero.
 * genotypes with many individuals are resampled using a Gaussian distribution
 */
int haploid_gt_dis::resample(double n)
{
	double pop_size;
	if (n<1.0) pop_size=population_size;
	else pop_size=n;

	population.set_state(HC_FUNC);
	double threshold_HG_CONTINUOUS=double(HG_CONTINUOUS)/pop_size, norm;
	norm=0;
	for (int i=0; i<(1<<number_of_loci); i++){
		if (population.func[i]<threshold_HG_CONTINUOUS)
		{
			population.func[i]=double(gsl_ran_poisson(rng, pop_size*population.func[i]))/pop_size;
		}
		else
		{
			population.func[i]+=double(gsl_ran_gaussian(rng, sqrt(population.func[i]/pop_size)));
		}
		norm+=population.func[i];
	}
	if (norm<HG_NOTHING){
		return HG_EXTINCT;
	}
	else population.scale(1.0/norm);
	return 0;
}

/*
 * calculate the distribution of mutants and update the population distribution
 */
int haploid_gt_dis::mutate()
{
	int locus;
	mutants.set_state(HC_FUNC);
	population.set_state(HC_FUNC);
	for (int i=0; i<(1<<number_of_loci); i++){
		mutants.func[i]=0;
		for (locus=0; locus<number_of_loci; locus++)
		{
			if (i&(1<<locus)){
				mutants.func[i]+=mutation_rates[0][locus]*population.func[i-(1<<locus)]-mutation_rates[1][locus]*population.func[i];
			}else{
				mutants.func[i]+=mutation_rates[1][locus]*population.func[i+(1<<locus)]-mutation_rates[0][locus]*population.func[i];
			}
		}
	}
	for (int i=0; i<(1<<number_of_loci); i++){
		population.func[i]+=mutants.func[i];
		//cout <<mutation_rate<<"  "<<mutants.func[i]<<"  "<<population.func[i]<<endl;
	}
	return 0;
}

/*
 * Calculate the distribution of recombinants and update the population, in case of free
 * recombinations, a fraction is replaced (outcrossing_rate), in case of general recombination
 * the entire population is replaced, i.e. obligate mating.
 */
int haploid_gt_dis::recombine(){
	int err;
	err=calculate_recombinants();
	population.set_state(HC_FUNC);
	if (free_recombination){
		for (int i=0; i<(1<<number_of_loci); i++){
			population.func[i]+=outcrossing_rate*(recombinants.func[i]-population.func[i]);
		}
	}else{
		for (int i=0; i<(1<<number_of_loci); i++){
			population.func[i]=recombinants.func[i];
		}
	}
	return err;
}

/*
 * determine how to calculate the recombinants and call the relevant routine
 */
int haploid_gt_dis::calculate_recombinants()
{
	if (free_recombination) return calculate_recombinants_free();
	else return calculate_recombinants_general();
}

/*
 * calculate the recombinant distribution for the free recombination case
 * almost the same as for the more general case below, but kept separate for
 * performance reasons - this is the most expensive part (3^L)
 */
int haploid_gt_dis::calculate_recombinants_free()
{
	int i,j,k, maternal_alleles, paternal_alleles, count;
	population.fft_func_to_coeff();
	recombinants.set_state(HC_COEFF);

	//normalization of the distribution
	recombinants.coeff[0]=1.0/(1<<number_of_loci);
	//loop of all coefficients of the distribution of recombinants
	for (i=1; i<(1<<number_of_loci); i++)
	{
		recombinants.coeff[i]=0;
		//loop over all possible partitions of the loci s1..sk in R^(k)_s1..sk to mother and father
		for (j=0; j<(1<<recombinants.order[i]); j++)
		{
			count=0;
			maternal_alleles=0;
			paternal_alleles=0;
			//build the integers to pull out the maternal and paternal moments
			for (k=0; k<number_of_loci; k++)
			{
				if (i&(1<<k))
				{
					if (j&(1<<count)) maternal_alleles+=(1<<k);
					else paternal_alleles+=(1<<k);
					count++;
				}
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
 * @brief calculate recombinants in the general case
 *
 * Calculate the distribution after recombination assumed in random mating with
 * pairs sampled with replacement.
 */
int haploid_gt_dis::calculate_recombinants_general()
{
	int i,j,k, maternal_alleles, paternal_alleles, count;
	population.fft_func_to_coeff();
	recombinants.set_state(HC_COEFF);

	//normalization of the distribution
	recombinants.coeff[0]=1.0/(1<<number_of_loci);
	//cout <<0<<"  "<<recombinants.coeff[0]<<endl;
	//loop of all coefficients of the distribution of recombinants
	for (i=1; i<(1<<number_of_loci); i++)
	{
		recombinants.coeff[i]=0;
		//loop over all possible partitions of the loci s1..sk in R^(k)_s1..sk to mother and father
		for (j=0; j<(1<<recombinants.order[i]); j++)
		{
			count=0;
			maternal_alleles=0;
			paternal_alleles=0;
			//build the integers to pull out the maternal and paternal moments
			for (k=0; k<number_of_loci; k++)
			{
				if (i&(1<<k))
				{
					if (j&(1<<count)) maternal_alleles+=(1<<k);
					else paternal_alleles+=(1<<k);
					count++;
				}
			}
			//add this particular contribution to the recombinant distribution
			recombinants.coeff[i]+=recombination_patters[i][j]*population.coeff[maternal_alleles]*population.coeff[paternal_alleles];
			//cout <<i<<"  "<<recombinants.coeff[i]<<"  "<<population.coeff[paternal_alleles]<<endl;
		}
		//normalize: the factor 1<<number_of_loci is due to a peculiarity of the the fft algorithm
		recombinants.coeff[i]*=(1<<(number_of_loci));
		//cout <<i<<"  "<<recombinants.coeff[i]<<endl;
	}
	//backtransform to genotype representation
	recombinants.fft_coeff_to_func();
	return 0;
}

/*
 * four routines that set the mutation rates, depending on what type of input is given
 */

/*
 * uniform mutation rate, input one double number
 */
int haploid_gt_dis::set_mutation_rate(double m){
	if (mem){
		for (int fb=0; fb<2; fb++){
			for (int locus=0; locus<number_of_loci; locus++){
				mutation_rates[fb][locus]=m;
			}
		}
		return 0;
	}else{
		cerr<<"haploid_gt_dis::set_mutation_rate(): allocate memory first!\n";
		return HG_MEMERR;
	}
}
/*
 * uniform mutation ratem backward and forward different, input two double number
 */
int haploid_gt_dis::set_mutation_rate(double mforward, double mbackward){
	if (mem){
		for (int locus=0; locus<number_of_loci; locus++){
			mutation_rates[0][locus]=mforward;
			mutation_rates[1][locus]=mbackward;
		}
		return 0;
	}else{
		cerr<<"haploid_gt_dis::set_mutation_rate(): allocate memory first!\n";
		return HG_MEMERR;
	}
}
/*
 * locus specific mutation rate, backward forward same, input one double array
 */
int haploid_gt_dis::set_mutation_rate(double* m){
	if (mem){
		for (int locus=0; locus<number_of_loci; locus++){
			mutation_rates[0][locus]=m[locus];
			mutation_rates[1][locus]=m[locus];
		}
		return 0;
	}else{
		cerr<<"haploid_gt_dis::set_mutation_rate(): allocate memory first!\n";
		return HG_MEMERR;
	}
}
/*
 * locus specific mutation rate, backward forward differebt, 2d double array
 */
int haploid_gt_dis::set_mutation_rate(double** m) {
	if (mem){
		for (int fb=0; fb<2; fb++){
			for (int locus=0; locus<number_of_loci; locus++){
				mutation_rates[fb][locus]=m[fb][locus];
			}
		}
		return 0;
	}else{
		cerr<<"haploid_gt_dis::set_mutation_rate(): allocate memory first!\n";
		return HG_MEMERR;
	}
}


/**
 * @brief calculate recombination patterns
 *
 * @param rec_rates a vector of recombination rates. The first entry should be large for linear chromosomes.
 *
 * @returns zero if successful, error codes otherwise (e.g. out of memory)
 *
 * A routine the calculates the probability of all possible recombination patters and
 * subpatterns thereof from a vector of recombination rates (rec_rates) passed as argument.
 * It allocated the memory (3^L) and calculates the entire distribution.
 *
 * The first entry is the recombination rate before the first locus, i.e. it should be large >50
 * for linear chromosomes. all other entries are recombination rates between successive loci.
 *
 */
int haploid_gt_dis::set_recombination_rates(double *rec_rates) {
	double err=0;
	int i, spin;
	//check whether the memory is already allocated, do so if not
	if (free_recombination==true)
	{
		int temp;
		int *nspins;	//temporary variables the track the number of ones in the binary representation of i
		nspins=new int [1<<number_of_loci];
		if (nspins==NULL) {
			cerr<<"haploid_gt_dis::set_recombination_rates(): Can not allocate memory!"<<endl;
			return HG_MEMERR;
		}
		spin=-1;
		nspins[0]=0;
		//allocate space for all possible subsets of loci
		recombination_patters=new double* [1<<number_of_loci];
		recombination_patters[0]=new double	[1];
		//loop over all possible locus subsets and allocate space for all
		//possible ways to assign the subset to father and mother (2^nspins)
		for (i=1; i<(1<<number_of_loci); i++){
			if (i==(1<<(spin+1))) spin++;
			temp=1+nspins[i-(1<<spin)];	//the order of coefficient k is 1+(the order of coefficient[k-2^spin])
			nspins[i]=temp;
			//all possible ways to assign the subset to father and mother (2^nspins)
			recombination_patters[i]=new double [(1<<temp)];
			if (recombination_patters[i]==NULL) err+=1;
		}
		delete [] nspins;
	}
	if (err==0){	//if memory allocation has been successful, calculate the probabilities of recombination
		int strand=0,newstrand, locus, set_size, i, subset, rec_pattern, marg_locus, higher_order_subset, higher_order_rec_pattern;
		double *rptemp;
		double sum=0;
		int strandswitches;
		//calculate the probabilities of different cross over realizations
		//the constrained of even number of crossovers is fulfilled automatically
		for (i=0; i<(1<<number_of_loci); i++){
			recombination_patters[(1<<number_of_loci)-1][i]=1.0;
			strand=(i&(1<<(number_of_loci-1)))>0?1:0;
			strandswitches=0;
			for (locus=0; locus<number_of_loci; locus++)
			{
				newstrand=((i&(1<<locus))>0)?1:0;
				if (strand==newstrand) recombination_patters[(1<<number_of_loci)-1][i]*=(0.5*(1.0+exp(-2.0*rec_rates[locus])));
				else {
					recombination_patters[(1<<number_of_loci)-1][i]*=(0.5*(1.0-exp(-2.0*rec_rates[locus])));
					strandswitches++;
				}
				strand=newstrand;
			}
			if (strandswitches%2) recombination_patters[(1<<number_of_loci)-1][i]=0;
			sum+=recombination_patters[(1<<number_of_loci)-1][i];
		}
		for (i=0; i<(1<<number_of_loci); i++){
			recombination_patters[(1<<number_of_loci)-1][i]/=sum;
		}
		//loop over set of spins of different size, starting with 11111101111 type patters
		//then 11101110111 type patterns etc. first loop is over different numbers of ones, i.e. spins
		for (set_size=number_of_loci-1; set_size>=0; set_size--)
		{
			//loop over all 2^L binary patterns
			for (subset=0; subset<(1<<number_of_loci); subset++)
			{
				//if correct number of ones... (its the same in every hypercube...)
				if (fitness.order[subset]==set_size)
				{
					marg_locus=-1; //determine the first zero, i.e. a locus that can be used to marginalize
					for (locus=0; locus<number_of_loci; locus++)
					{
						if ((subset&(1<<locus))==0)
							{marg_locus=locus; break;}
					}
					//a short hand for the higher order recombination pattern, from which we will marginalize
					higher_order_subset=subset+(1<<marg_locus);
					rptemp=recombination_patters[higher_order_subset];
					//loop over all pattern of the length set_size and marginalize
					//i.e. 111x01011=111001011+111101011
					for (rec_pattern=0; rec_pattern<(1<<set_size); rec_pattern++){
						higher_order_rec_pattern=(rec_pattern&((1<<marg_locus)-1))+((rec_pattern&((1<<set_size)-(1<<marg_locus)))<<1);
						recombination_patters[subset][rec_pattern]=rptemp[higher_order_rec_pattern]+rptemp[higher_order_rec_pattern+(1<<marg_locus)];
					}
				}
			}
		}
		free_recombination=false;
		return 0;
	}else{
		cerr <<"haploid_gt_dis::set_recombination_rates(): cannot allocate memory for recombination patterns!"<<endl;
		return HG_MEMERR;
	}
}

/*
 * calculate the genotype entropy and return
 */
double haploid_gt_dis::genotype_entropy(){
	double S=0;
	if (population.get_state()==HC_COEFF) population.fft_coeff_to_func();
	for (int i=0; i<(1<<number_of_loci); i++){
		S-=population.func[i]*log(population.func[i]);
	}
	return S;
}

/*
 * calculate the allele entropy and return
 * it has be made sure that the population.fft_func_to_coeff() was called
 */
double haploid_gt_dis::allele_entropy(){
	double SA=0;
	if (population.get_state()==HC_FUNC) population.fft_func_to_coeff();
	for (int locus=0; locus<number_of_loci; locus++){
		SA-=0.5*(1.0+population.coeff[(1<<locus)])*log(0.5*(1.0+population.coeff[(1<<locus)]));
		SA-=0.5*(1.0-population.coeff[(1<<locus)])*log(0.5*(1.0-population.coeff[(1<<locus)]));
	}
	return SA;
}

/*
 * calculate the fitness mean and variance and return
 * it has be made sure that the population.fft_func_to_coeff() was called
 */
stat_t haploid_gt_dis::get_fitness_statistics(){
	double mf=0, sq=0, temp;
	if (population.get_state()==HC_COEFF) population.fft_coeff_to_func();
	for (int locus=0; locus<1<<number_of_loci; locus++){
		temp=population.get_func(locus)*fitness.get_func(locus);
		mf+=temp;
		sq+=temp*temp;
	}
	return stat_t(mf, sq-mf);
}


/*
 * debugging routine: calculates the distribution of recombinants explicitly and
 * compares the result to the recombinant distribution obtained via fourier transform
 */
int haploid_gt_dis::test_recombinant_distribution(){
	double *test_rec;
	double dev=0;
	//allocate memory for the recombinant distribution calculated step-by-step
	test_rec=new double [(1<<number_of_loci)];
	int mother, father;
	//calculate recombinants the efficient way
	calculate_recombinants();
	//now calculate the recombinant distribution from pairs of parents.
	int gt1, gt2, rec_pattern;
	if (free_recombination){
		for (gt1=0; gt1<(1<<number_of_loci); gt1++){	//target genotype
			test_rec[gt1]=0.0;							//initialize
			//loop over all recombination patters (equal probability)
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
		for (gt1=0; gt1<(1<<number_of_loci); gt1++){
			test_rec[gt1]=0.0;
			for (rec_pattern=0; rec_pattern<(1<<number_of_loci); rec_pattern++){
				for (gt2=0; gt2<(1<<number_of_loci); gt2++){
					mother=(gt1&(rec_pattern))+(gt2&(~rec_pattern));
					father=(gt1&(~rec_pattern))+(gt2&(rec_pattern));
					//contribution is weighted by the probability of this particular recombination pattern
					//this got calculated and stored in recombination_patters[(1<<number_of_loci)-1]
					test_rec[gt1]+=recombination_patters[(1<<number_of_loci)-1][rec_pattern]*population.func[mother]*population.func[father];
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

/*
 * Debugging routine: produces random genotypes configurations and test whether they recombine right
 */
int haploid_gt_dis::test_recombination(double *rec_rates){

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



int haploid_gt_dis::mutation_drift_equilibrium(double **mu){
	set_mutation_rate(mu);
	//init population and recombination rates
	double *af=new double[number_of_loci];;
	double *recrates=new double[number_of_loci];;
	for (int i=0; i<number_of_loci; i++){
		af[i]=0;
		recrates[i]=10;
	}
	init_frequencies(af);
	//allocate histograms to store allele frequency distributions
	gsl_histogram **mutfreq=new gsl_histogram* [number_of_loci];
	for (int locus=0; locus<number_of_loci; locus++){
		mutfreq[locus]=gsl_histogram_alloc(100);
		gsl_histogram_set_ranges_uniform(mutfreq[locus], -1,1);
	}

	//equilibrate for 2N generations
	for (int gen=0; gen<2*population_size; gen++){
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
			theorynorm[locus]+=pow(0.5*(1+0.5*(upper+lower)), 2*population_size*mu[0][locus]-1)*pow(0.5*(1-0.5*(upper+lower)), 2*population_size*mu[1][locus]-1);
		}
	}
	for (int i=0; i<100; i++){
		gsl_histogram_get_range(mutfreq[0], i, &lower, &upper);
		cout <<setw(15)<<0.5*(upper+lower);
		for (int locus=0; locus<number_of_loci; locus++){
			cout <<setw(15)<<gsl_histogram_get(mutfreq[locus], i)/histogramnorm[locus]
					<<setw(15)<<pow(0.5*(1+0.5*(upper+lower)), 2*population_size*mu[0][locus]-1)*pow(0.5*(1-0.5*(upper+lower)), 2*population_size*mu[1][locus]-1)/theorynorm[locus];
		}
		cout <<endl;
	}
	return 0;
}
