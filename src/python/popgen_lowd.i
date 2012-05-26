/* renames and ignores */
%ignore hypercube_lowd;
%ignore haploid_lowd_test;

/* additional helper functions */
%pythoncode {
def binarify(gt, L=0):
    '''Transform an integer into a binary sequence on the L hypercube.'''
    import numpy as np
    if not L:
        L=1
        while gt > ((1<<L) - 1):
            L += 1
    return np.array(map(lambda l: bool(gt&(1<<(L-l-1))),range(L)))


def integerify(b):
    '''Transform a binary sequence on the HC into an integer.'''
    import numpy as np
    L = len(b)
    a = [(1<<(L-l-1)) for l in xrange(L)]
    return np.dot(b,a)
}

/**** HAPLOID_GT_DIS ****/
%define DOCSTRING_HAPLOID_GT_DIS
"Class for low-dimensional population genetics (short genomes ~20 loci).

This class is the main object for simulating the evolution of populations with
a few loci (less than ~20). The class offers a number of functions, but an
example will explain the basic idea:

#####################################
#   EXAMPLE SCRIPT                  #
#####################################
import numpy as np
import matplotlib.pyplot as plt
import PopGenLib as h

c = h.haploid_lowd(5)
c.set_genotypes([0, 2], [300, 700]) 
c.set_fitness_additive([0.02,0.03,0.04,0.02, -0.03])
c.evolve(10)
c.plot_diversity_histogram()
plt.show()
#####################################

An effective way to discover all available methods is to import PopGenLib from
an interactive shell (e.g. iPython), create a population as above, and use TAB
autocompletion:

In [1]: import PopGenLib as h
In [2]: c = h.haploid_lowd(5, 2000)
In [3]: c.      <--- TAB
"
%enddef
%feature("autodoc", DOCSTRING_HAPLOID_GT_DIS) haploid_lowd;
%extend haploid_lowd {

/* constructor */
%define DOCSTRING_HAPLOID_GT_DIS_INIT
"Construct a low-dimensional population with certain parameters.

Parameters:
- L : number of loci (at least 1)
- rng_seed : seed for the random number generator    
"
%enddef
%feature("autodoc", DOCSTRING_HAPLOID_GT_DIS_INIT) haploid_lowd;
%exception haploid_lowd {
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
        sprintf(buffer,"haploid_lowd: L = %d, N = %d", $self->L(), $self->N());
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"haploid_lowd(%d, %5.2e)", $self->L(), $self->N());
        return &buffer[0];
}

/* TODO: ignore hypercubes for now */
%ignore fitness;
%ignore population;


/* read only parameters */
%ignore L;
%ignore N;
%rename (_get_number_of_loci) get_number_of_loci;
%rename (_get_population_size) get_population_size;
%rename (_get_generation) get_generation;
%pythoncode {
L = property(_get_number_of_loci)
N = property(_get_population_size)
number_of_loci = property(_get_number_of_loci)
population_size = property(_get_population_size)
generation = property(_get_generation)
}


/* initialize frequencies */
%ignore set_allele_frequencies;
int _set_allele_frequencies(int DIM1, double *IN_ARRAY1, unsigned long N) {return $self->set_allele_frequencies(IN_ARRAY1, N);}
%pythoncode {
def set_allele_frequencies(self, frequencies, N=1000):
    '''Initialize the population in linkage equilibrium with allele frequencies.

    Parameters:
    - frequencies: an array of length L with all allele frequencies
    - N: the carrying capacity (target population size)

    *Note*: the latter parameter is only used for resampling and has therefore
    no crucial effect on the speed of the simulation.
    '''
    if len(frequencies) != self.L:
        raise ValueError('The input array of allele frequencies has the wrong length.')
    if self._set_allele_frequencies(frequencies, N):
        raise RuntimeError('Error in the C++ function.')
}

/* initialize genotypes */
%ignore set_genotypes;
%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* indices), (int len2, double* vals)};
int _set_genotypes(int len1, double* indices, int len2, double* vals) {
        vector<index_value_pair_t> gt;
        index_value_pair_t temp;
        for(size_t i = 0; i != len1; i++) {
                temp.index = (int)indices[i];
                temp.val = vals[i];
                gt.push_back(temp);
        }
        return $self->set_genotypes(gt);
}
%clear (int len1, double* indices);
%clear (int len2, double* vals);
%pythoncode {
def set_genotypes(self, indices, counts):
    '''Initialize population with fixed counts for specific genotypes.

    Parameters:
    - indices: list of genotypes to set (e.g. 0 --> 00...0, L-1 --> 11...1)
    - counts: list of counts for those genotypes

    *Note*: the population size and the carrying capacity are set as the sum of the counts.
    *Note*: you can use Python binary notation for the indices, e.g. 0b0110 is 6.
    '''
    import numpy as np
    indices = np.asarray(indices, float)
    counts = np.asarray(counts, float)
    if len(indices) != len(counts):
        raise ValueError('Indices and counts must have the same length')
    if self._set_genotypes(indices, counts):
        raise RuntimeError('Error in the C++ function.')
}

/* set recombination rates */
%typemap(in) double *rec_rates {
        /* Ensure input is a Python sequence */
        PyObject *tmplist = PySequence_Fast($input, "I expected a sequence");
        unsigned long L = PySequence_Length(tmplist);

        /* Get circular and L properties from the class (we are in the Python world here) */
        bool circular = (bool)PyInt_AsLong(PyObject_GetAttrString($self, "circular"));
        long Lint = PyInt_AsLong(PyObject_GetAttrString($self, "L"));

        /* Check lengths */
        if((!(circular)) && (L != Lint - 1)) {
                PyErr_SetString(PyExc_ValueError, "Expecting an array of length L-1.");
                SWIG_fail;
        }        
        if((circular) && (L != Lint)) {
                PyErr_SetString(PyExc_ValueError, "Expecting an array of length L.");
                SWIG_fail;
        } 

        /* Create C array from Python list */
        $1 = new double[L];
        double tmpdouble;
        for(size_t i=0; i < L; i++) {
                tmpdouble = (double)PyFloat_AsDouble(PySequence_Fast_GET_ITEM(tmplist, i));
                if (tmpdouble < 0) {
                        PyErr_SetString(PyExc_ValueError,"Expecting a sequence of positive floats");
                        SWIG_fail;
                }
                $1[i] = tmpdouble;
        }
}
%typemap(freearg) double *rec_rates {
  if($1) delete[] $1;
}

/* mutation rate(s) */
%rename (_get_mutation_rate) get_mutation_rate;
%pythoncode {
def get_mutation_rate(self, locus=None, direction=None):
    '''Get one or several mutation rates.

    Parameters:
    - locus: get only the mutation rate(s) of this locus
    - direction: get only the forward or backward mutation rate(s)

    Returns:
    - the mutation rate(s) requensted

    *Note*: if the mutation rates for all loci and/or directions are the same,
    this function will try to be smart and give you the answer you are looking for.
    In case of doubt, you will get a matrix (L x 2) with the full mutation rate
    landscape.
    '''

    import numpy as np
    if locus is not None:
        if not np.isscalar(locus):
            raise TypeError('Please select a *single* locus or no locus at all.')
        if direction is not None:
            return self._get_mutation_rate(locus, direction)
        else:
            mrs = tuple([self._get_mutation_rate(locus, d) for d in [0,1]])
            if mrs[0] == mrs[1]:
                return mrs[0]
            else:
                return mrs
    else:
        if direction is not None:
            mrs = np.array([self._get_mutation_rate(l, direction) for l in xrange(self.L)])
            if len(np.unique(mrs)) == 1:
                return mrs[0]
            else:
                return mrs
        else:
            mrs = np.array([[self._get_mutation_rate(l, d) for l in xrange(self.L)] for d in [0,1]])
            if len(np.unique(mrs)) == 1:
                return mrs[0,0]
            else:
                return mrs
}

%ignore set_mutation_rate;
int _set_mutation_rate(double *IN_ARRAY2, int DIM1, int DIM2) {
        double ** mrs = new double*[2];
        for(size_t i = 0; i < 2; i++)
                mrs[i] = &(IN_ARRAY2[DIM2 * i]);
        int result = $self->set_mutation_rate(mrs);
        delete[] mrs;
        return result;
}
%pythoncode{        
def set_mutation_rate(self, rates, rates_back=None):
    '''Set the mutation rate.

    Parameters:
    - rates: if a double, the mutation rate at any locus in both directions
         or, if rates_back is not None, only in the forward direction
         if a vector, the mutation rate is specified for each locus, the same
         in both directions or, if rates_back is not None, only in the
         forward direction
    - rates_back: mutation rate in the backward direction (global or
         locus-specific)
    '''

    import numpy as np
    L = self.L
    if np.isscalar(rates):
        if rates_back is None:
            ratesm = np.repeat(rates, L * 2).reshape(2,L)
        else:
            ratesm = np.vstack([np.repeat(rates, L), np.repeat(rates_back, L)])
    elif (np.rank(rates) != 1) or ((rates_back is not None) and (np.rank(rates_back) != 1)):
        raise ValueError('Please input one/two numbers or arrays.')
    else:
        if rates_back is None:
            ratesm = np.vstack([rates, rates])
        else:
            ratesm = np.vstack([rates, rates_back])

    if self._set_mutation_rate(ratesm):
        raise RuntimeError('Error in the C++ function.')
}

/* genotype frequencies */
%pythoncode {
def get_genotype_frequencies(self):
    '''Get the frequency of each genotype.'''
    import numpy as np
    return np.array([self.get_genotype_frequency(l) for l in xrange(1<<self.L)])
}

/* allele frequencies */
%pythoncode {
def get_allele_frequencies(self):
    '''Get all allele frequencies'''
    import numpy as np
    return np.array([self.get_allele_frequency(l) for l in xrange(self.L)])
}

/* ignore tests (they work by now) */
%ignore test_recombinant_distribution();
%ignore test_recombination(double *rec_rates);
%ignore mutation_drift_equilibrium(double** mutrates);

/* random sampling */
%pythoncode {
def random_clones(self, n_sample):
    '''Get random clones according to their frequencies.'''
    import numpy as np
    counts = np.random.multinomial(n_sample, self.get_genotype_frequencies())
    ind = counts.nonzero()[0]
    counts = counts[ind]
    sample = np.concatenate([np.repeat(ind[i], counts[i]) for i in xrange(len(ind))])
    np.random.shuffle(sample)
    return sample
}

/* get fitnesses of all individuals */
void _get_fitnesses(int DIM1, double* ARGOUT_ARRAY1) {
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = $self->get_fitness(i);
}
%pythoncode {
def get_fitnesses(self):
    '''Get the fitness of all possible genotypes.'''
    return self._get_fitnesses(1<<self.L)
}

/* divergence/diversity/fitness distributions and plot (full Python implementations) */
%pythoncode {
def get_fitness_histogram(self, n_sample=1000, **kwargs):
    '''Get the histogram of the fitness in the population.'''
    import numpy as np

    # Random sample
    gt = self.random_clones(n_sample)

    # Calculate fitness
    fit = np.array([self.get_fitness(gt[i]) for i in xrange(n_sample)])

    return np.histogram(fit, bins=bins, **kwargs)


def plot_fitness_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot the histogram of the fitness in the population.'''
    import numpy as np
    import matplotlib.pyplot as plt

    # Random sample
    gt = self.random_clones(n_sample)

    # Calculate fitness
    fit = np.array([self.get_fitness(gt[i]) for i in xrange(n_sample)])

    # Plot
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Fitness histogram')
        axis.set_xlabel('Fitness')
    axis.hist(fit, **kwargs)


def get_divergence_statistics(self, n_sample=1000):
    '''Get the mean and variance of the divergence in the population.'''
    import numpy as np
    L = self.L

    # Random sample
    gt = self.random_clones(n_sample)

    # Calculate divegence
    div = np.array([binarify(gt[i], L).sum() for i in xrange(n_sample)], int)

    return stat(div.mean(), div.var())


def get_divergence_histogram(self, bins=10, n_sample=1000, **kwargs):
    '''Get the histogram of the divergence in the population.'''
    import numpy as np
    L = self.L

    # Random sample
    gt = self.random_clones(n_sample)

    # Calculate divergence
    div = np.array([binarify(gt[i], L).sum() for i in xrange(n_sample)], int)

    return np.histogram(div, bins=bins, **kwargs)


def plot_divergence_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot the histogram of the divergence in the population.'''
    import numpy as np
    import matplotlib.pyplot as plt
    L = self.L

    # Random sample
    gt = self.random_clones(n_sample)

    # Calculate divegence
    div = np.array([binarify(gt[i], L).sum() for i in xrange(n_sample)], int)

    # Plot
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Divergence histogram')
        axis.set_xlabel('Divergence')
    
    if 'bins' not in kwargs:
        kwargs['bins'] = np.arange(10) * max(1, (div.max() + 1 - div.min()) / 10) + div.min()
    axis.hist(div, **kwargs)


def get_diversity_statistics(self, n_sample=1000):
    '''Get the mean and variance of the diversity in the population.'''
    import numpy as np
    L = self.L

    # Random sample
    gt1 = self.random_clones(n_sample)
    gt2 = self.random_clones(n_sample)

    # Calculate diversity
    div = np.array([binarify(gt1[i] ^ gt2[i], L).sum() for i in xrange(n_sample)], int)

    return stat(div.mean(), div.var())


def get_diversity_histogram(self, bins=10, n_sample=1000, **kwargs):
    '''Get the histogram of the diversity in the population.'''
    import numpy as np
    L = self.L

    # Random sample
    gt1 = self.random_clones(n_sample)
    gt2 = self.random_clones(n_sample)

    # Calculate diversity
    div = np.array([binarify(gt1[i] ^ gt2[i], L).sum() for i in xrange(n_sample)], int)

    # Calculate histogram
    return np.histogram(div, bins=bins, **kwargs)


def plot_diversity_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot the histogram of the diversity in the population.'''
    import numpy as np
    import matplotlib.pyplot as plt
    L = self.L

    # Random sample
    gt1 = self.random_clones(n_sample)
    gt2 = self.random_clones(n_sample)

    # Calculate diversity
    div = np.array([binarify(gt1[i] ^ gt2[i], L).sum() for i in xrange(n_sample)], int)

    # Plot
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Diversity histogram')
        axis.set_xlabel('Diversity')
    
    if 'bins' not in kwargs:
        kwargs['bins'] = np.arange(10) * max(1, (div.max() + 1 - div.min()) / 10) + div.min()
    axis.hist(div, **kwargs)
}

/* set fitness landscape */
%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* indices), (int len2, double* vals)};
int _set_fitness_func(int len1, double* indices, int len2, double* vals) {
        vector<index_value_pair_t> iv;
        index_value_pair_t temp;
        for(size_t i = 0; i != len1; i++) {
                temp.index = (int)indices[i];
                temp.val = vals[i];
                iv.push_back(temp);
        }
        return ($self->fitness).init_list(iv);
}
%clear (int len1, double* indices);
%clear (int len2, double* vals);
%pythoncode {
def set_fitness_function(self, indices, vals):
    '''Set the fitness landscape at single points.

    Parameters:
    - indices: genotype to which the fitness values will be assigned
    - vals: fitness values to assign
    '''
    import numpy as np
    indices = np.asarray(indices, float)
    vals = np.asarray(vals, float)
    if len(indices) != len(vals):
        raise ValueError('Indices and values must have the same length')
    if self._set_fitness_func(indices, vals):
        raise RuntimeError('Error in the C++ function.')
}

/* set additive fitness component */
%exception set_fitness_additive {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
void set_fitness_additive(int DIM1, double* IN_ARRAY1) {
        if(DIM1 != $self->L())
                PyErr_Format(PyExc_ValueError, "The array had a wrong length.");
        if (($self->fitness).additive(IN_ARRAY1))
                PyErr_Format(PyExc_RuntimeError, "Error in the C++ function.");
}
} /* extend haploid_lowd */
