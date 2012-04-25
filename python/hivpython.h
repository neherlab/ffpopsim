/**
 * @file hivpopulation.h
 * @brief Header file for a typical HIV population (subclass of haploid_clone)
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-23
 */
#ifndef HIVPYTHON_H_
#define HIVPYTHON_H_

#include <string>
#include <sstream>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#define HIVPOP_VERBOSE 0
#define HIVPOP_BADARG -1354341
#define NOTHING 1e-10
#define HIVGENOME 10000


/**
 * @brief HIV population with facultative drug treatment
 *
 * This class exemplifies the haploid_clone base class. It mainly adds one trait,
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
 *
 * This class is an interface class for the Python bindings (like a dressing screen).
 */
class hivpython {
private:
	//random number generator
	double treatment;
//	gsl_rng* rng;
//	int seed;
//	using haploid_clone::set_up;	// only the new set_up functino is allowed, lest people mess with the genome length
//public:
	// constructors/destructors
	hivpython();
	virtual ~hivpython();
//	int set_up(int N_in, int rng_seed=0, double mutrate=3e-5, double coinfection_rate=1e-2, double crossover_rate=1e-3);
//
//	// treatment (set/get)
//	void set_treatment(double t){treatment=t; update_traits(); update_fitness();}
//	double get_treatment() {return treatment;}
//	void calc_fitness_from_traits(clone_t *tempgt) {tempgt->fitness = tempgt->trait[0] + treatment * tempgt->trait[1];}
//
//	// stream I/O
//	int read_selection_coefficients(istream &model);
//	int read_resistance_coefficients(istream &model);
//	int write_genotypes(ostream &out, int sample_size, string gt_label, int start=0, int length=0);

};

#endif /* HIVPYTHON_H_ */
