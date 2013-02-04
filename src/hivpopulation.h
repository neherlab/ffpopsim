/**
 * @file hivpopulation.h
 * @brief Header file for a typical HIV population (subclass of haploid_highd)
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-04-23
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
#ifndef HIVPOPULATION_H_
#define HIVPOPULATION_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include "ffpopsim_highd.h"

#define HIVPOP_VERBOSE 0
#define HIVPOP_BADARG -1354341
#define NOTHING 1e-10
#define HIVGENOME 10000

// HIV genes
#define ENV_START 7000
#define ENV_END 8000


/**
 * @brief HIV gene.
 *
 * Attributes:
 * - start: the starting position of the gene
 * - stop: the (last position + 1) of the gene
 */
struct hivgene {
	unsigned int start;
	unsigned int end;
	hivgene(unsigned int start_in=0, unsigned int end_in=HIVGENOME);
};

/**
 * @brief HIV population with facultative drug treatment
 *
 * This class exemplifies the haploid_highd base class. It mainly adds one trait,
 * "treatment", which is the same for all individuals and represents the presence
 * or absence of drug treatment (in a continuous manner, \f$0 \leq \f$ treatment
 * \f$\leq 1\f$).
 *
 * The replication capacity in absence of drug is encoded in the first trait. The
 * drug resistance phenotype is represented by the second trait. Fitness is
 * computed from traits as follows:
 *
 * f[trait] = trait[0] + treatment * trait[1]
 *
 * Moreover, this class fixes the length of the genome to exactly 10000 sites.
 */
class hivpopulation : public haploid_highd {
public:
	// constructors/destructors
	hivpopulation(int N=0, int rng_seed=0, double mutation_rate=3e-5, double coinfection_rate=1e-2, double crossover_rate=1e-3);
	virtual ~hivpopulation();

	// genes
	hivgene env;

	// treatment (set/get)
	void set_treatment(double t){treatment=t; update_traits(); update_fitness();}
	double get_treatment() {return treatment;}

	// stream I/O
	int read_replication_coefficients(istream &model);
	int read_resistance_coefficients(istream &model);
	int write_genotypes(ostream &out_genotypes, int sample_size, string gt_label="", int start=0, int length=0);

protected:
	// fitness landscape
	virtual void calc_individual_fitness_from_traits(clone_t *tempgt);

private:
	//random number generator
	double treatment;
	gsl_rng* rng;
	int seed;

};

#endif /* HIVPOPULATION_H_ */
