/*
 *  hypercube_function.cpp
 *  Haploids
 *
 *  Created by Richard Neher on 9/19/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "popgen.h"

hypercube_function::hypercube_function()
{
mem=false;
if (HC_VERBOSE) cerr<<"hypercube_function::hypercube_function(): constructing...!\n";
}

//set up and allocate memory, perform basic consistency checks on input.
int hypercube_function::set_up(int dim_in, int s)
{
	if (dim_in>0)
	{
		if (HC_VERBOSE) cerr<<"hypercube_function::set_up(): setting up...!";
		dim=dim_in;
		mem=false;
		hcube_allocated=false;
		epistatic_std=0;
		hypercube_mean=0.0;
		coefficients_epistasis.clear();
		coefficients_single_locus.clear();
		if (HC_VERBOSE) cerr<<"done.\n";
		if (s==0) s=time(NULL);
		seed=s;
		return allocate_mem();
	}
	else
	{
		cerr <<"hypercube_function: need positive dimension!\n";
		return HC_BADARG;
	}
}

//destructor
hypercube_function::~hypercube_function()
{
	if (HC_VERBOSE) cerr<<"hypercube_function::~hypercube_function(): destructing...!\n";
	if (mem) free_mem();
}

//allocate the necessary memory
int hypercube_function::allocate_mem()
{
	if (HC_VERBOSE) cerr<<"hypercube_function::allocate_mem(): allocating memory...";
	if (mem)
	{
		cerr <<"hypercube_function::allocate_mem(): memory already allocated, freeing and reallocating ...!\n";
		free_mem();
	}

	//a random number generator used for the coefficients and the like
	rng=gsl_rng_alloc(RNG);
	gsl_rng_set(rng, seed);
	mem=true;
	if (HC_VERBOSE) cerr<<"done.\n";
	return 0;
}

//free the memory
int hypercube_function::free_mem()
{
	if (!mem)
	{
		cerr <<"hypercube_function::free_mem(): no memory allocated...!\n";
		return 0;
	}
 	mem=false;
 	gsl_rng_free(rng);
 	return 0;
}


double hypercube_function::get_func(boost::dynamic_bitset<> *gt)
{
	if (HC_VERBOSE) cerr<<"fluct_hypercube::get_func()\n";
	double result=hypercube_mean;
	int sign=1, locus;
	//first order contributions
	for (unsigned int c=0; c<coefficients_single_locus.size();c++){
		if ((*gt)[coefficients_single_locus[c].locus]) result+=coefficients_single_locus[c].value;
		else result-=coefficients_single_locus[c].value;
	}
	//interaction contributions
	//if (HC_VERBOSE)
	for (unsigned int c=0; c<coefficients_epistasis.size();c++){
		sign=1;
		for (locus=0; locus<coefficients_epistasis[c].order; locus++){
			if (!(*gt)[coefficients_epistasis[c].loci[locus]]) sign*=-1;
		}
		result+=sign*coefficients_epistasis[c].value;
	}

	//calculate the random fitness part
	//TODO random epistasis is disabled for now since the old partition in to integers no
	//longer exists. This should be easy to remedy in a similar way the random reassortment
	//patterns are done in haploid_clone
	int gt_seed=0;
	/*if (epistatic_std>0)
	{
		//calculate the seed for the random number generator
		for (int i=0; i<n_o_w; i++) gt_seed+=(*gt)[i];
		//add a gaussion random number to the fitness from the rng seeded with the genoytpe
		gsl_rng_set(rng,gt_seed+rng_offset);
		result+=gsl_ran_gaussian(rng,epistatic_std);
	}*/
	return result;
}

/*
 * old version of the previous function which takes an array of integers rather than a bitset
 * the second argument specifies how many bits of each integer are used to store the genotype
 */
double hypercube_function::get_func_words(int *gt, int n_o_w)
{
	if (HC_VERBOSE) cerr<<"fluct_hypercube::get_func_words()\n";
	double result=hypercube_mean;
	int sign=1, locus;
	//first order contributions
	for (unsigned int c=0; c<coefficients_single_locus.size();c++){
		if ((gt[coefficients_single_locus[c].word]&(1<<coefficients_single_locus[c].bit))) result+=coefficients_single_locus[c].value;
		else result-=coefficients_single_locus[c].value;
	}
	//interaction contributions
	for (unsigned int c=0; c<coefficients_epistasis.size();c++){
		sign=1;
		for (locus=0; locus<coefficients_epistasis[c].order; locus++){
			if (!(gt[coefficients_epistasis[c].words[locus]]&(1<<coefficients_epistasis[c].bits[locus]))) sign*=-1;
		}
		result+=sign*coefficients_epistasis[c].value;
	}

	//calculate the random fitness part
	int gt_seed=0;
	if (epistatic_std>0)
	{
		//calculate the seed for the random number generator
		for (int i=0; i<n_o_w; i++) gt_seed+=gt[i];
		//add a gaussion random number to the fitness from the rng seeded with the genoytpe
		gsl_rng_set(rng,gt_seed+rng_offset);
		result+=gsl_ran_gaussian(rng,epistatic_std);
	}
	return result;
}

/*
 * function that takes a double and assigns it to a trait coefficients with index set
 * specified in vector loci. It is overloaded, see below.
 * this version assumes that genotypes are stored in integer arrays and take the size of that
 * array and the number of bits in each integer as additional arguments
 */
int hypercube_function::add_coefficient(double value, vector <int> loci, int n_o_w, int* n_o_b)
{
	if (loci.size()>1){
		coeff temp_coeff(value, loci);
		coefficients_epistasis.push_back(temp_coeff);
		int temp_locus=0,w;
		for (int locus=0; locus<coefficients_epistasis.back().order; locus++){
			temp_locus=coefficients_epistasis.back().loci[locus];
			if (temp_locus!=loci[locus]) {cerr <<"hypercube_function::add_coefficient: should never happen!"; return HC_BADARG;}
			w=0;
			while (w<n_o_w and temp_locus>n_o_b[w]) {temp_locus-=n_o_b[w]; w++;}
			coefficients_epistasis.back().words[locus]=w;
			coefficients_epistasis.back().bits[locus]=temp_locus;
			if (w==n_o_w) {cerr <<"hypercube_function::add_coefficient: should never happen!"; return HC_BADARG;}
		}
	}else if (loci.size()==1){
		coeff_single_locus temp_coeff(value, loci[0]);
		coefficients_single_locus.push_back(temp_coeff);
		int temp_locus=0,w;
		temp_locus=coefficients_single_locus.back().locus;
		if (temp_locus!=loci[0]) {cerr <<"hypercube_function::add_coefficient: should never happen!"; return HC_BADARG;}
		w=0;
		while (w<n_o_w and temp_locus>=n_o_b[w]) {temp_locus-=n_o_b[w]; w++;}
		coefficients_single_locus.back().word=w;
		coefficients_single_locus.back().bit=temp_locus;
		//cerr <<"hypercube_function::add_coefficient: locus"<<coefficients_single_locus.back().locus<<" word "<<coefficients_single_locus.back().word<<" bit "<<coefficients_single_locus.back().bit<<endl;
		if (w==n_o_w) {cerr <<"hypercube_function::add_coefficient: should never happen!"; return HC_BADARG;}
	}else{ hypercube_mean=value;}
	return 0;
}

/*
 * version of the above for bitsets. the word and bits attribute of the coefficient are not set
 */
int hypercube_function::add_coefficient(double value, vector <int> loci)
{
	if (loci.size()>1){
		coeff temp_coeff(value, loci);
		coefficients_epistasis.push_back(temp_coeff);
	}else if (loci.size()==1){
		coeff_single_locus temp_coeff(value, loci[0]);
		coefficients_single_locus.push_back(temp_coeff);
	}else{ hypercube_mean=value;}
	return 0;
}

/*
 * reset the value of an additive coefficient, requires the index of the coefficient
 * and the locus it is supposed to refer to
 */
int hypercube_function::set_additive_coefficient(double value, int lindex, int expected_locus)
{
	//check whether the locus is what was expected
	if (coefficients_single_locus[lindex].locus==expected_locus){
		coefficients_single_locus[lindex].value=value;
		return 0;
	}else if (expected_locus>-1){
		cerr <<"hypercube_function::set_additive_coefficient: coefficient[locus] does not match locus: "<<coefficients_single_locus[lindex].locus<<"  "<<expected_locus<<endl;
		return HC_BADARG;
	}
}

/*
 * set the standard deviation of random epistatic fitness, as of now, only works for genotype
 * specified by words rather than bitsets
 */
//TODO update comment
int hypercube_function::set_random_epistasis_strength(double sigma)
{
	if (sigma>0){epistatic_std=sigma; return 0;}
	else{
		cerr<<"hypercube_function::set_random_epistasis_strength(): Epistasis strength has to be positive!"<<endl;
		return HC_BADARG;
	}
}
/*
 * loop over the additive coefficients and return the value once the locus is found.
 */
double hypercube_function::get_additive_coefficient(int locus){
	int l=0;
	while(l<coefficients_single_locus.size()){
		if (coefficients_single_locus[l].locus==locus) return coefficients_single_locus[l].value;
		l++;
	}
	return 0.0;
}

