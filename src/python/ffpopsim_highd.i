/* renames and ignores */
%ignore coeff_t;
%ignore coeff_single_locus_t;
%ignore hypercube_highd;

/**** CLONE_T ****/
%rename (clone) clone_t;
%extend clone_t {

/* string representations */
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"clone: %d traits, genome size = %d", ($self->trait).size(), ($self->genotype).size());
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"clone");
        return &buffer[0];
}


/* traits */
%rename (_trait) trait;
int _get_number_of_traits() {
        return ($self->trait).size();
}

void _get_trait(int DIM1, double* ARGOUT_ARRAY1) {
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = ($self->trait)[i];
}
%pythoncode {
number_of_traits = property(_get_number_of_traits)
@property
def trait(self):
        return self._get_trait(self.number_of_traits)
}

/* genotype */
%typemap(in) boost::dynamic_bitset<> genotype_in (boost::dynamic_bitset<> temp) {
        /* Ensure input is a Python sequence */
        PyObject *tmplist = PySequence_Fast($input, "I expected a sequence");
        unsigned long L = PySequence_Length(tmplist);

        /* Create boost::dynamic_bitset from Python list */
        temp.resize(L);
        long tmplong;
        for(size_t i=0; i < L; i++) {
                tmplong = PyInt_AsLong(PySequence_Fast_GET_ITEM(tmplist, i));
                if(tmplong < 0) {
                        PyErr_SetString(PyExc_ValueError, "Expecting an array of bool.");
                        SWIG_fail;
                }
                temp[i] = (bool)tmplong; 
        }      
        $1 = temp;
}

%ignore genotype;
int _get_genotype_length() {return ($self->genotype).size();}
void _get_genotype(int DIM1, short* ARGOUT_ARRAY1) {
        for(size_t i=0; i < ($self->genotype).size(); i++) ARGOUT_ARRAY1[i] = ($self->genotype)[i];
}

void _set_genotype(boost::dynamic_bitset<> genotype_in) {$self->genotype = genotype_in;}

%pythoncode {
@property
def genotype(self):
        import numpy as np
        return np.array(self._get_genotype(self._get_genotype_length()), bool)


@genotype.setter
def genotype(self, genotype):
        self._set_genotype(genotype)
}
} /* extend clone_t */

/**** HAPLOID_CLONE ****/
%define DOCSTRING_HAPLOID_CLONE
"Class for high-dimensional population genetics (genomes larger than ~20 loci).

This class is the main object for simulating the evolution of populations with
many loci (more than ~20). The class offers a number of functions, but an
example will explain the basic idea:

#####################################
#   EXAMPLE SCRIPT                  #
#####################################
import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h

c = h.haploid_highd(5000, 2000)
c.set_genotypes() 
c.mutation_rate = 0.01
c.evolve(10)
c.plot_divergence_histogram()
plt.show()
#####################################

An effective way to discover all available methods is to import FFPopSim from
an interactive shell (e.g. iPython), create a population as above, and use TAB
autocompletion:

In [1]: import FFPopSim as h
In [2]: c = h.haploid_highd(5000, 2000)
In [3]: c.      <--- TAB


Populations can have a number of phenotypic traits that concur to the fitness
of each individual. The function that calculates fitness from the phenotype is
(TODO: should be) described by the user, as well as all the full genotype-
phenotype map.

*Note*: fitness is not a phenotypic trait directly, but rather a function of _all_
phenotypic traits together. However, in case 
"
%enddef
%feature("autodoc", DOCSTRING_HAPLOID_CLONE) haploid_highd;

%extend haploid_highd {

/* constructor */
%define DOCSTRING_HAPLOID_CLONE_INIT
"Construct a high-dimensional population with certain parameters.

Parameters:
- L     length of the genome(number of loci)
- rng_seed      seed for the random generator. If zero (default) pick a random number
- number_of_traits      number of phenotypic traits
"
%enddef
%feature("autodoc", DOCSTRING_HAPLOID_CLONE_INIT) haploid_highd;

/* ignore problematic stuff */
%ignore traits;
/* TODO: implement a Python list as a surrogate for a pointer, and expose hypercube_highd */

/* constructor */
%exception haploid_highd {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Construction impossible. Please check input args.");
                SWIG_fail;
        }
}

/* string representations */
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"haploid_highd: L = %d, N = %d", $self->L(), $self->N());
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"haploid_highd(%d, %5.2e)", $self->L(), $self->N());
        return &buffer[0];
}

/* read only parameters */
%ignore L;
%ignore N;
%rename (_get_number_of_loci) get_number_of_loci;
%rename (_get_population_size) get_population_size;
%rename (_get_generation) get_generation;
%rename (_get_number_of_clones) get_number_of_clones;
%rename (_get_number_of_traits) get_number_of_traits;
%rename (max_fitness) get_max_fitness;
%pythoncode {
L = property(_get_number_of_loci)
N = property(_get_population_size)
number_of_loci = property(_get_number_of_loci)
population_size = property(_get_population_size)
generation = property(_get_generation)
number_of_clones = property(_get_number_of_clones)
number_of_traits = property(_get_number_of_traits)
}

/* initalize frequencies */
%ignore set_allele_frequencies;
int _set_allele_frequencies(double *IN_ARRAY1, int DIM1, int n_o_genotypes) {return $self->set_allele_frequencies(IN_ARRAY1, n_o_genotypes);}
%pythoncode {
def set_allele_frequencies(self, frequencies, N=1000):
    '''Initialize the population according to the given allele frequencies.

    Parameters:
    - frequencies: an array of length L with all allele frequencies
    - N: the carrying capacity (target population size)
    '''
    if len(frequencies) != self.L:
            raise ValueError('Please input an L dimensional list of allele frequencies.')
    if self._set_allele_frequencies(frequencies, N):
        raise RuntimeError('Error in the C++ function.')
}

/* initialize genotypes */
%ignore set_genotypes;
%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* genotypes), (int len2, double* vals)};
int _set_genotypes(int len1, double* genotypes, int len2, double* vals) {
        /* We use a flattened array */
        len1 /= len2;
        vector<genotype_value_pair_t> gt;
        genotype_value_pair_t temp;
        for(size_t i = 0; i != len1; i++) {
                temp.genotype = boost::dynamic_bitset<>(len1);
                for(size_t j=0; j < len1; j++)
                        temp.genotype[j] = (bool)genotypes[i * len1 + j];
                temp.val = vals[i];
                gt.push_back(temp);
        }
        return $self->set_genotypes(gt);
}
%clear (int len1, double* genotypes);
%clear (int len2, double* vals);
%pythoncode {
def set_genotypes(self, genotypes, counts):
    '''Initialize population with fixed counts for specific genotypes.

    Parameters:
    - indices: list of genotypes to set (e.g. 0 --> 00...0, L-1 --> 11...1)
    - counts: list of counts for those genotypes

    *Note*: the population size and the carrying capacity are set as the sum of the counts.
    *Note*: you can use Python binary notation for the indices, e.g. 0b0110 is 6.
    '''
    import numpy as np
    genotypes = np.array(genotypes, float, copy=False, ndmin=2)
    counts = np.asarray(counts, float)
    if len(genotypes) != len(counts):
        raise ValueError('Indices and counts must have the same length')
    if self._set_genotypes(genotypes.flatten(), counts):
        raise RuntimeError('Error in the C++ function.')
}



/* evolve */
%rename (_evolve) evolve;
%pythoncode{
def evolve(self, gen=1):
        '''Evolve for some generations.

        Parameters:
        - gen: number of generations
        '''

        if self._evolve(gen):
                raise RuntimeError('Error in the C++ function.')
        else:
                self.calc_stat()
}


/* get allele frequencies */
void _get_allele_frequencies(double* ARGOUT_ARRAY1, int DIM1) {
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = $self->get_allele_frequency(i);
}
%pythoncode {
def get_allele_frequencies(self): return self._get_allele_frequencies(self.L)
}

/* get genotypes */
void get_genotype(unsigned int i, short* ARGOUT_ARRAY1, int DIM1) {
        boost::dynamic_bitset<> newgt = (*($self->current_pop))[i].genotype;
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = newgt[i];
}
%pythoncode {
def get_genotypes(self, ind=None):
        '''Get genotypes of the population.

        Parameters:
        - ind: if a scalar, a single genotype corresponding to clone ind is returned.
               otherwise, several genotypes are returned (default: all)
        '''
        import numpy as np
        L = self.number_of_loci
        if np.isscalar(ind):
                return np.array(self.get_genotype(ind, L), bool)

        if ind is None:
                ind = xrange(self.number_of_clones)
        genotypes = np.zeros((len(ind), L), bool)
        for i, indi in enumerate(ind):
                genotypes[i] = self.get_genotype(indi, L)
        return genotypes
}

/* set trait/fitness coefficients */
%exception clear_fitness {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Fitness depends only on traits, not on the genome directly.");
                SWIG_fail;
        }
}

/* create trait (fitness) landscapes */
void clear_trait(unsigned int traitnumber) {
        /* check trait number */
        if(traitnumber >= $self->get_number_of_traits())
                throw HIVPOP_BADARG;
        ($self->trait)[traitnumber].reset();
}

/* get single locus effects */
void _get_additive_trait(double* ARGOUT_ARRAY1, int DIM1, int traitnumber=0) {
        /* check trait number */
        if(traitnumber >= $self->get_number_of_traits())
                throw HIVPOP_BADARG;

        /* Initialize to zero */
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = 0;

        /* Add any coefficient you found */
        hypercube_highd *trait = &(($self->trait)[traitnumber]);
        coeff_single_locus_t * coeff;
        for(size_t i=0; i < trait->coefficients_single_locus.size(); i++) {
                coeff = &(trait->coefficients_single_locus[i]);
                ARGOUT_ARRAY1[coeff->locus] += coeff->value;
        }
}
%pythoncode{
def get_additive_trait(self, traitnumber=0):
        '''Get an array with the additive part of a trait for all loci.'''
        return self._get_additive_trait(self.L, traitnumber)
}

/* set single locus effects */
void set_additive_trait(int DIM1, double* IN_ARRAY1, int traitnumber=0) {
        /* check trait number */
        if(traitnumber >= $self->get_number_of_traits())
                throw HIVPOP_BADARG;
        /* check length of vector */
        if(DIM1 != $self->L())
                throw HIVPOP_BADARG;

        vector <int> loci;
        for(size_t i = 0; i < DIM1; i++) {
                loci.push_back(i);
                $self->add_trait_coefficient(IN_ARRAY1[i], loci, traitnumber);
                loci.clear();
        }
}


/* add_trait_coefficient: this is implemented as a conversion from a Python vector of loci to std::vector */
%typemap(in) vector<int> loci (std::vector<int> temp) {
        /* Ensure input is a Python sequence */
        PyObject *tmplist = PySequence_Fast($input, "I expected a sequence");
        unsigned long L = PySequence_Length(tmplist);

        /* Create std::vector from Python list */
        temp.reserve(L);
        long tmplong;
        for(size_t i=0; i < L; i++) {
                tmplong = PyInt_AsLong(PySequence_Fast_GET_ITEM(tmplist, i));
                if(tmplong < 0) {
                        PyErr_SetString(PyExc_ValueError, "Expecting an array of positive integers (the loci).");
                        SWIG_fail;
                }
                temp.push_back((int)tmplong); 
        }      
        $1 = temp;
}


/* get fitnesses of all clones */
void _get_fitnesses(int DIM1, double* ARGOUT_ARRAY1) {
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = $self->get_fitness(i);
}
%pythoncode {
def get_fitnesses(self):
        '''Get the fitness of all clones.'''
        return self._get_fitnesses(self.number_of_clones)
}

/* Hamming distance (full Python reimplementation) */
%ignore distance_Hamming;
%pythoncode {
def distance_Hamming(self, clone_gt1, clone_gt2, chunks=None, every=1):
        import numpy as np
        if np.isscalar(clone_gt1):
                genotypes = self.get_genotypes((clone_gt1, clone_gt2))
                clone_gt1 = genotypes[0]
                clone_gt2 = genotypes[1]

        if chunks is not None:
                ind = np.zeros(clones.shape[1], bool)
                for chunk in chunks:
                        inde = np.arange(chunk[1] - chunk[0])
                        inde = inde[(inde % every) == 0] + chunk[0]
                        ind[inde] = True
                clone_gt1 = clone_gt1[ind]
                clone_gt2 = clone_gt2[ind]
        return (clone_gt1 != clone_gt2).sum()
}

/* get random clones/genotypes */
%ignore random_clones(unsigned int n_o_individuals, vector <int> *sample);
%pythoncode {
def random_genomes(self, n):
        import numpy as np
        L = self.number_of_loci
        genotypes = np.zeros((n, L), bool)
        for i in xrange(genotypes.shape[0]):
                genotypes[i] = self.get_genotype(self.random_clone(), L)
        return genotypes
}
void random_clones(int DIM1, unsigned int * ARGOUT_ARRAY1) {
        vector <int> sample = vector <int>(0);
        int err = $self->random_clones(DIM1, &sample);
        if(!err)
                for(size_t i=0; i < DIM1; i++)
                        ARGOUT_ARRAY1[i] = sample[i];
}

/* divergence/diversity/fitness distributions and plot */
%ignore get_divergence_histogram;
%ignore get_diversity_histogram;
%ignore get_fitness_histogram;
%pythoncode {
def get_fitness_histogram(self, bins=10, n_sample=1000, **kwargs):
    '''Calculate the fitness histogram.
    Parameters:
    - bins: number or array of bins to be used in the histogram (see also numpy.histogram)
    - n_sample: number of individuals to sample
    '''
    import numpy as np
	#self.produce_random_sample(n_sample)
    fit = [self.get_fitness(self.random_clone()) for i in xrange(n_sample)]
    h = np.histogram(fit, bins=bins, **kwargs)
    return h
    
    
def plot_fitness_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot a distribution of fitness in the population.
    Parameters:
    - axis: an axis to use. A new figure is created by default
    - n_sample: number of individuals to sample
    
    Note: any other keyword arguments (bins, color, etc.) are passed verbatim to the
          matplotlib.pyplot.hist function
    '''
    import matplotlib.pyplot as plt
	#self.produce_random_sample(n_sample)
    fit = [self.get_fitness(self.random_clone()) for i in xrange(n_sample)]
    
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Fitness histogram')
        axis.set_xlabel('Fitness')
    return axis.hist(fit, **kwargs)
    
    
def get_divergence_histogram(self, bins=10, chunks=None, every=1, n_sample=1000, **kwargs):
    '''Get the divergence histogram restricted to those chunks of the genome.

    Parameters:
    - bins: number or array of bins to be used in the histogram (see also numpy.histogram)
    - chunks: restrict analysis to some chunk in the genome. It must be an n x 2 matrix with
              the initial and (final+1) positions of the chunks
    - every: restrict analysis to every X positions. For instance, if every third site is neutral,
             this argument can be used to only look at those neutral sites
    - n_sample: number of individuals to sample
    '''
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
    return np.histogram(div, bins=bins, **kwargs)
    
    
def plot_divergence_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot the divergence histogram.

    Parameters:
    - axis: an axis to use. A new figure is created by default
    - n_sample: number of individuals to sample
    
    Note: any other keyword arguments (bins, color, etc.) are passed verbatim to the
          matplotlib.pyplot.hist function
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    genotypes = self.random_genomes(n_sample)
    div = genotypes.sum(axis=1)
     
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Divergence histogram')
        axis.set_xlabel('Divergence')
    
    if 'bins' not in kwargs:
        kwargs['bins'] = np.arange(10) * max(1, (div.max() + 1 - div.min()) / 10) + div.min()
    return axis.hist(div, **kwargs)
    
    
def get_diversity_histogram(self, bins=10, chunks=None, every=1, n_sample=1000, **kwargs):
    '''Get the diversity histogram restricted to those chunks of the genome.

    Parameters:
    - bins: number or array of bins to be used in the histogram (see also numpy.histogram)
    - chunks: restrict analysis to some chunk in the genome. It must be an n x 2 matrix with
              the initial and (final+1) positions of the chunks
    - every: restrict analysis to every X positions. For instance, if every third site is neutral,
             this argument can be used to only look at those neutral sites
    - n_sample: number of individuals to sample
    '''
    import numpy as np
    
    # Check chunks
    if chunks is not None:
        chunks = np.asarray(chunks)
        if (np.rank(chunks) != 2) or (chunks.shape[1] != 2):
            raise ValueError('Please input an N x 2 matrix with the chunks initial and (final+1) positions')
    
    # Get the random genotypes
    genotypes = self.random_genomes(2 * n_sample)
    
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
    return np.histogram(div, bins=bins, **kwargs)


def plot_diversity_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot the diversity histogram.

    Parameters:
    - axis: an axis to use. A new figure is created by default
    - n_sample: number of individuals to sample
    
    Note: any other keyword arguments (bins, color, etc.) are passed verbatim to the
          matplotlib.pyplot.hist function
    '''
    import matplotlib.pyplot as plt
    import numpy as np
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
    return axis.hist(div, **kwargs)
}

} /* extend haploid_highd */
