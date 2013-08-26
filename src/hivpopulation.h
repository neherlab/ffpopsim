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
//  - coordinates refer to HXB2
//  - HXB2 is not exactly 10000 bases long, but this makes no relevant difference
#define GAG_START 789
#define GAG_END 2292
#define POL_START 2087
#define POL_END 5096
#define ENV_START 6314
#define ENV_END 8795
#define NEF_START 8796
#define NEF_END 9417
#define VIF_START 5040
#define VIF_END 5619
#define VPR_START 5558
#define VPR_END 5850
#define VPU_START 6061
#define VPU_END 6310
#define REV1_START 5969
#define REV1_END 6045
#define REV2_START 8378
#define REV2_END 8653
#define TAT1_START 5830
#define TAT1_END 6045
#define TAT2_START 8378
#define TAT2_END 8469



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
	unsigned int second_start;
	unsigned int second_end;
	hivgene(unsigned int start_in=0, unsigned int end_in=HIVGENOME,
		unsigned int second_start_in=0, unsigned int second_end_in=0);
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
	hivgene gag;
	hivgene pol;
	hivgene env;
	hivgene nef;
	hivgene vif;
	hivgene vpu;
	hivgene vpr;
	hivgene tat;
	hivgene rev;

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
