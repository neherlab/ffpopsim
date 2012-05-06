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

/* Numpy magic to output arrays */
/*%include "pyfragments.swg"*/
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

    // random clone (combine this with get_genotype to get the genotype)
    int random_clone();
    
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
def random_genomes(self, n):
    import numpy as np
    genotypes = np.zeros((n, self.get_number_of_loci()), bool)
    for i in xrange(genotypes.shape[0]):
        genotypes[i] = self.get_genotype(self.random_clone())
    return genotypes


def get_fitnesses(self):
    return self._get_fitnesses(self.get_number_of_clones())


def get_fitness_histogram(self, bins=10, n_sample=1000, **kwargs):
    '''Calculate the fitness histogram.'''
    import numpy as np
    fit = [self.get_fitness(self.random_clone()) for i in xrange(n_sample)]
    h = np.histogram(fit, bins=bins, **kwargs)
    return h


def plot_fitness_histogram(self, axis=None, **kwargs):
    '''Plot a distribution of fitness in the population.'''
    import matplotlib.pyplot as plt
    fit = self.get_fitnesses();

    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Fitness histogram')
        axis.set_xlabel('Fitness')

    axis.hist(fit, **kwargs)


def get_divergence_histogram(self, bins=10, chunks=None, every=1, n_sample=1000, **kwargs):
    '''Get the divergence histogram restricted to those chunks of the genome.'''
    import numpy as np

    # Check chunks
    if chunks is not None:
        chunks = np.asarray(chunks)
        if (np.rank(chunks) != 2) or (chunks.shape[1] != 2):
            raise ValueError('Please input an N x 2 matrix with the chunks initial and (final+1) positions')

    # Get the random genotypes
    genotypes = self.random_genomes(n_sample)

    # Restrict to the chunks
    if chunks is not None:
        ind = np.zeros(genotypes.shape[1], bool)
        for chunk in chunks:
            inde = np.arange(chunk[1] - chunk[0])
            inde = inde[(inde % every) == 0] + chunk[0]
            ind[inde] = True
        genotypes = genotypes[:,ind]

    # Calculate divergence
    div = genotypes.sum(axis=1)

    # Calculate histogram
    h = np.histogram(div, bins=bins, **kwargs)
    return h 


def plot_divergence_histogram(self, axis=None, **kwargs):
    import matplotlib.pyplot as plt
    import numpy as np
    n_sample = 1000
    genotypes = self.random_genomes(n_sample)
    div = genotypes.sum(axis=1)
 
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Divergence histogram')
        axis.set_xlabel('Divergence')

    if 'bins' not in kwargs:
        kwargs['bins'] = np.arange(10) * max(1, (div.max() + 1 - div.min()) / 10) + div.min()

    axis.hist(div, **kwargs)


def get_diversity_histogram(self, bins=10, chunks=None, every=1, n_sample=1000, **kwargs):
    '''Get the diversity histogram restricted to those chunks of the genome.'''
    import numpy as np

    # Check chunks
    if chunks is not None:
        chunks = np.asarray(chunks)
        if (np.rank(chunks) != 2) or (chunks.shape[1] != 2):
            raise ValueError('Please input an N x 2 matrix with the chunks initial and (final+1) positions')

    # Get the random genotypes
    genotypes2 = self.random_genomes(2 * n_sample)
    for i in xrange(genotypes.shape[0]):
        genotypes[i] = self.get_genotype(self.random_clone())

    # Restrict to the chunks
    if chunks is not None:
        ind = np.zeros(genotypes.shape[1], bool)
        for chunk in chunks:
            inde = np.arange(chunk[1] - chunk[0])
            inde = inde[(inde % every) == 0] + chunk[0]
            ind[inde] = True
        genotypes = genotypes[:,ind]

    # Calculate diversity
    genotypes1 = genotypes[:genotypes.shape[0] / 2]
    genotypes2 = genotypes[-genotypes1.shape[0]:]
    div = (genotypes1 != genotypes2).sum(axis=1)

    # Calculate histogram
    h = np.histogram(div, bins=bins, **kwargs)
    return h 


def plot_diversity_histogram(self, axis=None, **kwargs):
    import matplotlib.pyplot as plt
    import numpy as np
    n_sample = 1000
    genotypes1 = self.random_genomes(n_sample)
    genotypes2 = self.random_genomes(n_sample)
    div = (genotypes1 != genotypes2).sum(axis=1)

    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Diversity histogram')
        axis.set_xlabel('Diversity')

    if 'bins' not in kwargs:
        kwargs['bins'] = np.arange(10) * max(1, (div.max() + 1 - div.min()) / 10) + div.min()

    axis.hist(div, **kwargs)


def write_genotypes_compressed(self, filename, sample_size, gt_label='', start=0, length=0):
    '''Write genotypes into a compressed archive.'''
    import numpy as np 
    if length <= 0:
        length = self.get_number_of_loci() - start
    d = {}
    for i in xrange(sample_size):
        rcl = self.random_clone()
        d['>'+str(i)+'_GT-'+gt_label+'_'+str(rcl)] = self.get_genotype(rcl)[start:start+length]
    np.savez_compressed(filename, **d)    
}
}

%include "hivpython.h"

