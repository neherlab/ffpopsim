/**
 * @file hivpopulation.h
 * @brief Header file for a typical HIV population (subclass of haploid_highd)
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-04-23
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
	virtual void calc_individual_fitness_from_traits(clone_t *tempgt) {tempgt->fitness = tempgt->trait[0] + treatment * tempgt->trait[1];}

private:
	//random number generator
	double treatment;
	gsl_rng* rng;
	int seed;

};

#endif /* HIVPOPULATION_H_ */
