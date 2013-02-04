/*
 *  hypercube_highd.cpp
 *  Haploids
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
#include "ffpopsim_highd.h"

hypercube_highd::hypercube_highd()
{
mem=false;
if (HCF_VERBOSE) cerr<<"hypercube_highd::hypercube_highd(): constructing...!\n";
}

//constructor
hypercube_highd::hypercube_highd(int dim_in, int s) {
	if (HCF_VERBOSE) cerr<<"hypercube_highd::hypercube_highd(): constructing...!\n";
	set_up(dim_in, s);
}

//set up and allocate memory, perform basic consistency checks on input.
int hypercube_highd::set_up(int dim_in, int s) {
	if (dim_in>0) {
		if (HCF_VERBOSE) cerr<<"hypercube_highd::set_up(): setting up...!";
                // FIXME: we should really go for a global rng seeder for the whole library
		if (s==0) s = time(NULL);
		seed = s;
		if (HCF_VERBOSE) cerr<<"done.\n";

                dim = dim_in;
                mem = false;
                hcube_allocated = false;
                hypercube_mean = 0;
                epistatic_std = 0;

		return allocate_mem();
	} else {
		cerr <<"hypercube_highd: need positive dimension!\n";
		return HCF_BADARG;
	}
}

//destructor
hypercube_highd::~hypercube_highd() {
	if (HCF_VERBOSE) cerr<<"hypercube_highd::~hypercube_highd(): destructing...!\n";
	if (mem) free_mem();
}

//allocate the necessary memory
int hypercube_highd::allocate_mem() {
	if (HCF_VERBOSE) cerr<<"hypercube_highd::allocate_mem(): allocating memory...";
	if (mem) {
		cerr <<"hypercube_highd::allocate_mem(): memory already allocated, freeing and reallocating ...!\n";
		free_mem();
	}

	//a random number generator used for the coefficients and the like
	rng=gsl_rng_alloc(RNG);
	gsl_rng_set(rng, seed);
	rng_offset = gsl_rng_uniform_int(rng, 1000000);
        if (HCF_VERBOSE) {
        	cerr <<"hypercube_highd() random number seed: "<<seed<<endl;
	        cerr <<"hypercube_highd() random number offset: "<<rng_offset<<endl;
        }

	// static single locus coefficients
	coefficients_single_locus_static = vector<double>(dim, 0);

	mem=true;
	if (HCF_VERBOSE) cerr<<"done.\n";
	return 0;
}

//free the memory
int hypercube_highd::free_mem() {
	if (!mem) {
		cerr <<"hypercube_highd::free_mem(): no memory allocated...!\n";
		return 0;
	}
 	mem=false;
 	gsl_rng_free(rng);
 	return 0;
}

/**
 * @brief Get single value on the hypercube
 *
 * @param genotype Point of the hypercube
 *
 * @returns the value corresponding to that point
 */
double hypercube_highd::get_func(boost::dynamic_bitset<>& genotype) {
	if (HCF_VERBOSE) cerr<<"fluct_hypercube::get_func()"<<endl;
	double result=hypercube_mean;
	int sign, locus;
	// first order contributions
	for (coefficients_single_locus_iter = coefficients_single_locus.begin();
	     coefficients_single_locus_iter != coefficients_single_locus.end();
	     coefficients_single_locus_iter++) {
		if ((genotype)[coefficients_single_locus_iter->locus]) result+=coefficients_single_locus_iter->value;
		else result -= coefficients_single_locus_iter->value;
	}
	// interaction contributions
	for (coefficients_epistasis_iter = coefficients_epistasis.begin();
	     coefficients_epistasis_iter != coefficients_epistasis.end();
	     coefficients_epistasis_iter++) {
		sign=1;
		for (locus=0; locus < coefficients_epistasis_iter->order; locus++) {
			if (!(genotype[coefficients_epistasis_iter->loci[locus]])) sign*= -1;
		}
		result += sign * coefficients_epistasis_iter->value;
	}
	// calculate the random fitness part
	if (epistatic_std > HP_NOTHING) {
		int gt_seed=0;
		int word=0, locus, ii;
		// calculate the seed for the random number generator
		for (locus=0;locus<dim;) {
			word=0;
			ii=0;
			while (ii<WORDLENGTH and locus<dim) {
				if ((genotype)[locus]) word+=(1<<ii);
				ii++; locus++;
			}
			gt_seed+=word;
		}
		// add a gaussion random number to the fitness from the rng seeded with the genoytpe
		gsl_rng_set(rng,gt_seed+rng_offset);
		result+=gsl_ran_gaussian(rng,epistatic_std);
	}
	if (HCF_VERBOSE) cerr<<"...done"<<endl;
	return result;
}

/**
 * @brief Calculate difference between two hypercube points efficiently
 *
 * @param genotype1 first point on the hypercube
 * @param genotype2 second point on the hypercube
 * @param diffpos vector of positions at which they differ
 *
 * @returns the difference between the values, f(gt1) - f(gt2)
 */
double hypercube_highd::get_func_diff(boost::dynamic_bitset<>& genotype1, boost::dynamic_bitset<>& genotype2, vector<int> &diffpos) {
	if (HCF_VERBOSE) cerr<<"fluct_hypercube::get_func_diff()"<<endl;
	double result = 0;
	int locus;

	//in case of random epistasis, evaluating the difference does not make a lot of sense.
	if (epistatic_std>HP_NOTHING){
		double val1, val2;
		val1 = get_func(genotype1);
		val2 = get_func(genotype2);
		return val1-val2;
	}else{
		// first order contributions
		for(size_t i=0; i != diffpos.size(); i++) {
			locus = diffpos[i];
			if (genotype1[locus] and !genotype2[locus]) result += 2 * get_additive_coefficient(locus);
			else if (!genotype1[locus] and genotype2[locus]) result -= 2 * get_additive_coefficient(locus);
			else{ cerr<<"fluct_hypercube::get_func_diff(): Difference vector is screwed up"<<endl;}
		}
		// TODO: calculate epistasis more efficiently!
		// interaction contributions
		int sign1, sign2;
		for (coefficients_epistasis_iter = coefficients_epistasis.begin();
			 coefficients_epistasis_iter != coefficients_epistasis.end();
			 coefficients_epistasis_iter++) {
			sign1 = sign2 = 1;
			for (locus=0; locus < coefficients_epistasis_iter->order; locus++) {
				if (!(genotype1[coefficients_epistasis_iter->loci[locus]])) sign1 *= -1;
				if (!(genotype2[coefficients_epistasis_iter->loci[locus]])) sign2 *= -1;
			}
			if(sign1 != sign2)
				result += (sign1 - sign2) * coefficients_epistasis_iter->value;
		}

		if (HCF_VERBOSE) cerr<<"...done"<<endl;
		return result;
	}
}


/**
 * @brief: get the trait coefficient of a locus.
 *
 * Loop over the additive coefficients and return the value once the locus is found.
 */
double hypercube_highd::get_additive_coefficient(int locus){
//	int l=0;
//	while(l<coefficients_single_locus.size()) {
//		if (coefficients_single_locus[l].locus==locus)
//			return coefficients_single_locus[l].value;
//		l++;
//	}
//	return 0.0;
	return coefficients_single_locus_static[locus];
}


/**
 * @brief Reset the hypercube
 */
void hypercube_highd::reset() {
	reset_additive();
	epistatic_std = 0;
	coefficients_epistasis.clear();
}


/**
 * @brief Reset the additive part of the hypercube
 */
void hypercube_highd::reset_additive() {
	hypercube_mean = 0;
	coefficients_single_locus.clear();
	fill(coefficients_single_locus_static.begin(), coefficients_single_locus_static.end(), 0);
}


/**
 * @brief assign single- or multi-locus trait coefficients.
 *
 * It takes a double and assigns it to a trait coefficients with index set
 * specified in vector loci.
 */
int hypercube_highd::add_coefficient(double value, vector <int> loci)
{
	if (loci.size()>1) {
		coeff_t temp_coeff(value, loci);
		coefficients_epistasis.push_back(temp_coeff);
	} else if (loci.size()==1) {
		coeff_single_locus_t temp_coeff(value, loci[0]);
		coefficients_single_locus.push_back(temp_coeff);
                coefficients_single_locus_static[loci[0]] = value;
	} else {
		hypercube_mean=value;
	}
	return 0;
}

/**
 * @brief reset the value of an additive coefficient
 *
 * Note: it requires the index of the coefficient and the locus it is supposed to refer to.
 */
int hypercube_highd::set_additive_coefficient(double value, int lindex, int expected_locus)
// TODO: do we need this function at all?
{
	//check whether the locus is what was expected
	if (coefficients_single_locus_static[expected_locus]==0 and value==0){
		return 0;
	}
	else if (coefficients_single_locus[lindex].locus==expected_locus){
		coefficients_single_locus[lindex].value=value;
		coefficients_single_locus_static[expected_locus]=value;
		return 0;
	}else if (expected_locus>-1){
		cerr <<"hypercube_highd::set_additive_coefficient: coefficient[locus] does not match locus: "<<coefficients_single_locus[lindex].locus<<"  "<<expected_locus<<endl;
		return HCF_BADARG;
	}
	return 0;
}

/**
 * @brief set the standard deviation of random epistatic fitness
 */
int hypercube_highd::set_random_epistasis_strength(double sigma)
{
	if (sigma>0){epistatic_std=sigma; return 0;}
	else{
		cerr<<"hypercube_highd::set_random_epistasis_strength(): Epistasis strength has to be positive!"<<endl;
		return HCF_BADARG;
	}
}
