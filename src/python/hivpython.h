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
 * @brief HIV population with facultative drug treatment (Python2 bindings)
 *
 * This class is an interface class for the Python bindings (like a dressing screen).
 * It contains an hivpopulation instance as a private attribute, so that problematic
 * bindings are avoided, and exposes only a few, Python-friendly methods.
 *
 */
class hivpython {
private:
	// hivpopulation class that does the job
	hivpopulation p;

public:
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
