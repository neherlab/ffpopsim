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
 */
class hivpython: private hivpopulation {
private:

public:
	// constructors/destructors
	hivpython();
	~hivpython();

	// we mask some baseclass methods by equivalent ones that are described in the .cpp file
	// this way we can selectively include code into SWIG
	int set_up(int N_in, int rng_seed=0, double mutrate=3e-5, double coinfection_rate=1e-2, double crossover_rate=1e-3){return hivpopulation::set_up(N_in, rng_seed, mutrate, coinfection_rate, crossover_rate);}

	// evolve
	int evolve(int gen=1);

	// random genotypes
	void random_clone(unsigned short ARGOUT_ARRAY1[HIVGENOME]);
	// Note: one cannot draw several clone in one shot because of limitations in
	// the SWIG Numpy interface (might be possible in the future).

	// allele frequencies
	double get_allele_frequency(int l){return haploid_clone::get_allele_frequency(l);}
	void get_allele_frequencies(double ARGOUT_ARRAY1[HIVGENOME]);

	// treatment (set/get)
	void set_treatment(double t){hivpopulation::set_treatment(t);}
	double get_treatment(){return hivpopulation::get_treatment();}
	void calc_fitness_from_traits(clone_t *tempgt){hivpopulation::calc_fitness_from_traits(tempgt);}

	// stream I/O
	int read_selection_coefficients(char *model);
	int read_resistance_coefficients(char *model);
//	int write_genotypes(ostream &out, int sample_size, string gt_label, int start=0, int length=0);

};

#endif /* HIVPYTHON_H_ */
