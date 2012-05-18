/**
 * @file hivpopulation.h
 * @brief Header file for a typical HIV population (subclass of haploid_highd)
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-23
 */
#ifndef HIVPOPULATION_H_
#define HIVPOPULATION_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include "popgen_highd.h"

#define HIVPOP_VERBOSE 0
#define HIVPOP_BADARG -1354341
#define NOTHING 1e-10
#define HIVGENOME 10000

// HIV genes
#define ENV_START 7000
#define ENV_END 8001


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
	hivpopulation(int N_in=0, int rng_seed=0, double mutation_rate_in=3e-5, double coinfection_rate_in=1e-2, double crossover_rate_in=1e-3);
	virtual ~hivpopulation();
	int set_up(int N_in, int rng_seed=0, double mutation_rate_in=3e-5, double coinfection_rate_in=1e-2, double crossover_rate_in=1e-3);

	// genes
	hivgene env;

	// treatment (set/get)
	void set_treatment(double t){treatment=t; update_traits(); update_fitness();}
	double get_treatment() {return treatment;}

	// stream I/O
	int read_replication_coefficients(istream &model);
	int read_resistance_coefficients(istream &model);
	int write_genotypes(ostream &out, int sample_size, string gt_label, int start=0, int length=0);

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
