/**
* @file hivpopulation.i
* @brief Python 2 bindings of the hivpopulation class.
* @author Richard Neher, Boris Shraiman, Fabio Zanini
* @version 
* @date 2012-04-24
*/
%module(docstring="Model for an HIV population.") hivpython
/* Include in the wrap code */
%{
#define SWIG_FILE_WITH_INIT
#include "hivpython.h"
%}

%ignore check_input;

/* Numpy magic to output arrays */
%include "numpy.i"
%init %{
import_array();
%}

/*
%apply (int DIM1, double* INPLACE_ARRAY1) {(int gts2, double* frequencies)};
%apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int gts, int loci, double* hoprates)};
%apply (int DIM1, double* IN_ARRAY1) {(int loci, double* fitness_additive)};
*/


/**************************************************************
 * CODE TO BE WRAPPED
 **************************************************************
 * Note: only _declarations_ are needed, no implementation is required at all for
 * wrapping.
 */

/**** STAT_T ****/
struct stat_t {
    double mean;
    double variance;
};

/**** HYPERCUBE_FUNCTION ****/
/*TODO: add C++ wrappers for the pointers (set/get methods) */
class hypercube_function {
public:
    int dim;
    double hypercube_mean;
    vector <coeff_single_locus_t> coefficients_single_locus;
    vector <coeff_t> coefficients_epistasis;
    double epistatic_std;
    
    // setting up
    hypercube_function();
    hypercube_function(int dim_in, int s=0);
    virtual ~hypercube_function();
    int set_up(int dim_in,  int s=0);
    
    // methods
    unsigned int get_seed() {return seed;};
    double get_func(boost::dynamic_bitset<> *genotype);
    double get_additive_coefficient(int locus);
    int set_additive_coefficient(double value, int locus, int expected_locus=-1);
    int add_coefficient(double value, vector <int> loci);
    int set_random_epistasis_strength(double sigma);
};


/**** CLONE_T ****/
%rename (_trait) clone_t::trait;
%rename (_genotype) clone_t::genotype;
struct clone_t {
    boost::dynamic_bitset<> genotype;
    vector <double> trait;
    double fitness;
    int clone_size;
    clone_t(int n_traits);
};
/* C++ HELPER CODE */
%extend clone_t {
        int get_number_of_traits() {
                return ($self->trait).size();
        }
        void _get_trait(int DIM1, double* ARGOUT_ARRAY1) {
                for(size_t i=0; i<($self->trait).size(); i++)
                        ARGOUT_ARRAY1[i] = ($self->trait)[i];
        }
        void _get_genotype(unsigned short ARGOUT_ARRAY1[HIVGENOME]) {
            for(size_t i=0; i < ($self->genotype).size(); i++) ARGOUT_ARRAY1[i] = ($self->genotype)[i];
        }

/* PYTHON HELPER CODE */
        %pythoncode {
        @property
        def trait(self):
            return self._get_trait(self.get_number_of_traits())
        
        @property
        def genotype(self):
                return self._get_genotype()
        }
}


/**** HAPLOID_CLONE ****/
class haploid_clone {
public:
    // population parameters (read only)
    int get_generation();
    int get_number_of_loci();
    int get_pop_size();
    int get_number_of_clones();
    
    //evolve
    int bottleneck(int size_of_bottleneck);
    
    // population parameters (read/write)
    int target_pop_size;
    double mutation_rate;
    double outcrossing_probability;
    double crossover_rate;
    int recombination_model;
    bool circular;
    
    // allele frequencies
    double get_allele_frequency(int l);
    double get_pair_frequency(int locus1, int locus2);
    
    // fitness/phenotype readout
    double get_fitness(int n);
    double get_trait(int n, int t=0);
};


/**************************************************************
 * INCLUDE HEADER OF THE SUBCLASS 'AS IS'
 *************************************************************/
/* PYTHON HELPER CODE */
%rename (_get_fitnesses) hivpython::get_fitnesses(int DIM1, double* ARGOUT_ARRAY1);
%extend hivpython {
%pythoncode {
def get_fitnesses(self):
    return self._get_fitnesses(self.get_number_of_clones())
}
}

%include "hivpython.h"

