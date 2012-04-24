/*
 *  hypercube_function.cpp
 *  Haploids
 *
 *  Created by Richard Neher on 9/19/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "popgen.h"
#include "popgen_highd.h"

hypercube_function::hypercube_function()
{
mem=false;
if (HCF_VERBOSE) cerr<<"hypercube_function::hypercube_function(): constructing...!\n";
}

//constructor
hypercube_function::hypercube_function(int dim_in, int s)
{
	if (HCF_VERBOSE) cerr<<"hypercube_function::hypercube_function(): constructing...!\n";
	set_up(dim_in, s);
}

//set up and allocate memory, perform basic consistency checks on input.
int hypercube_function::set_up(int dim_in, int s)
{
	if (dim_in>0)
	{
		if (HCF_VERBOSE) cerr<<"hypercube_function::set_up(): setting up...!";
		dim=dim_in;
		mem=false;
		hcube_allocated=false;
		epistatic_std=0;
		hypercube_mean=0.0;
		coefficients_epistasis.clear();
		coefficients_single_locus.clear();
		if (HCF_VERBOSE) cerr<<"done.\n";
		if (s==0) s=time(NULL);
		seed=s;
		return allocate_mem();
	}
	else
	{
		cerr <<"hypercube_function: need positive dimension!\n";
		return HCF_BADARG;
	}
}

//destructor
hypercube_function::~hypercube_function()
{
	if (HCF_VERBOSE) cerr<<"hypercube_function::~hypercube_function(): destructing...!\n";
	if (mem) free_mem();
}

//allocate the necessary memory
int hypercube_function::allocate_mem()
{
	if (HCF_VERBOSE) cerr<<"hypercube_function::allocate_mem(): allocating memory...";
	if (mem)
	{
		cerr <<"hypercube_function::allocate_mem(): memory already allocated, freeing and reallocating ...!\n";
		free_mem();
	}

	//a random number generator used for the coefficients and the like
	rng=gsl_rng_alloc(RNG);
	gsl_rng_set(rng, seed);
	mem=true;
	if (HCF_VERBOSE) cerr<<"done.\n";
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


double hypercube_function::get_func(boost::dynamic_bitset<> *genotype)
{
	if (HCF_VERBOSE) cerr<<"fluct_hypercube::get_func()\n";
	double result=hypercube_mean;
	int sign=1, locus;
	//first order contributions
	for (unsigned int c=0; c<coefficients_single_locus.size();c++){
		if ((*genotype)[coefficients_single_locus[c].locus]) result+=coefficients_single_locus[c].value;
		else result-=coefficients_single_locus[c].value;
	}
	//interaction contributions
	//if (HCF_VERBOSE)
	for (unsigned int c=0; c<coefficients_epistasis.size();c++){
		sign=1;
		for (locus=0; locus<coefficients_epistasis[c].order; locus++){
			if (!(*genotype)[coefficients_epistasis[c].loci[locus]]) sign*=-1;
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
		for (int i=0; i<n_o_w; i++) gt_seed+=(*genotype)[i];
		//add a gaussion random number to the fitness from the rng seeded with the genoytpe
		gsl_rng_set(rng,gt_seed+rng_offset);
		result+=gsl_ran_gaussian(rng,epistatic_std);
	}*/
	return result;
}

/**
 * @brief: get the trait coefficient of a locus.
 *
 * Loop over the additive coefficients and return the value once the locus is found.
 */
double hypercube_function::get_additive_coefficient(int locus){
	int l=0;
	while(l<coefficients_single_locus.size()){
		if (coefficients_single_locus[l].locus==locus)
			return coefficients_single_locus[l].value;
		l++;
	}
	return 0.0;
}

/**
 * @brief assign single- or multi-locus trait coefficients.
 *
 * It takes a double and assigns it to a trait coefficients with index set
 * specified in vector loci.
 */
int hypercube_function::add_coefficient(double value, vector <int> loci)
{
	if (loci.size()>1) {
		coeff_t temp_coeff(value, loci);
		coefficients_epistasis.push_back(temp_coeff);
	} else if (loci.size()==1) {
		coeff_single_locus_t temp_coeff(value, loci[0]);
		coefficients_single_locus.push_back(temp_coeff);
	} else { hypercube_mean=value;
	}
	return 0;
}

/**
 * @brief reset the value of an additive coefficient
 *
 * Note: it requires the index of the coefficient and the locus it is supposed to refer to.
 */
int hypercube_function::set_additive_coefficient(double value, int lindex, int expected_locus)
{
	//check whether the locus is what was expected
	if (coefficients_single_locus[lindex].locus==expected_locus){
		coefficients_single_locus[lindex].value=value;
		return 0;
	}else if (expected_locus>-1){
		cerr <<"hypercube_function::set_additive_coefficient: coefficient[locus] does not match locus: "<<coefficients_single_locus[lindex].locus<<"  "<<expected_locus<<endl;
		return HCF_BADARG;
	}
}

/**
 * @brief set the standard deviation of random epistatic fitness
 */
int hypercube_function::set_random_epistasis_strength(double sigma)
{
	if (sigma>0){epistatic_std=sigma; return 0;}
	else{
		cerr<<"hypercube_function::set_random_epistasis_strength(): Epistasis strength has to be positive!"<<endl;
		return HCF_BADARG;
	}
}
