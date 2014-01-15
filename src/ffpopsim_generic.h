/**
 * @file popgen.h
 * @brief Header file with the classes and types provided with the library. 
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2010-10-27
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
#ifndef FFPOPGEN_GENERIC_H_
#define FFPOPGEN_GENERIC_H_

#include <time.h>
#include <cmath>
#include <vector>
#include <bitset>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>

#define MIN(a,b) (a<b)?a:b
#define MAX(a,b) (a>b)?a:b
#define RNG gsl_rng_taus2		//choose the random number generator algorithm, see http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html

#define FREE_RECOMBINATION 1
#define CROSSOVERS 2
#define SINGLE_CROSSOVER 3

using namespace std;

/**
 * @brief Pairs of an index and a value
 */
struct index_value_pair_t {
	size_t index;
	double val;
	index_value_pair_t(int index_in=0, double val_in=0) : index(index_in), val(val_in) {};
};

/**
 * @brief Pairs of a genotype and a value
 */
struct genotype_value_pair_t {
	boost::dynamic_bitset<> genotype;
	double val;
	genotype_value_pair_t(boost::dynamic_bitset<> genotype_in=boost::dynamic_bitset<>(0), double val_in=0) : genotype(genotype_in), val(val_in) {};
};

/**
 * @brief Structure for short summary statistics.
 */
struct stat_t {
	double mean;
	double variance;
	stat_t(double mean_in=0, double variance_in=0) : mean(mean_in), variance(variance_in) {};
};

#define SAMPLE_ERROR -12312154

/**
 * @brief Sample of any scalar property.
 * 
 * This class is used to store samples of scalar quantities used in the evolution of the population,
 * for instance fitness or allele frequencies. I enables simple manipulations (mean, variance, etc.).
 */
class sample {
public:
	int number_of_values;
	double *values;
	double mean;
	double variance;

	int bins;
	bool mem_dis;
	bool mem_values;
	gsl_histogram *distribution;

	bool with_range;
	double range_min;
	double range_max;

	sample();
	virtual ~sample();
	int set_up(int n);
	int set_distribution(int bins=100);
	void set_range(double min, double max) {range_min=min; range_max=max; with_range=true;}
	int calc_mean();
	int calc_variance();
	int calc_distribution();
	int print_distribution(ostream &out);
};

#endif /* FFPOPGEN_GENERIC_H_ */
