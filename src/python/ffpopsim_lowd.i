/**
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

/* renames and ignores */
%ignore hypercube_lowd;
%ignore haploid_lowd_test;

/*****************************************************************************/
/* additional helper functions                                               */
/*****************************************************************************/
%pythoncode
%{
def binarify(gt, L=0):
    '''Transform an integer into a binary sequence on the L hypercube.

    Parameters:
       - gt: integer representing a genotype
       - L: number of dimensions of the hypercube

    Returns:
       - genotype: bool vector representing the same genotype

    **Examples**:

    .. sourcecode:: ipython

       In [1]: binarify(3, 5)
       Out[1]: array([False, False, False,  True,  True], dtype=bool)

       In [2]: FFPopSim.binarify(0b11, 5)
       Out[2]: array([False, False, False,  True,  True], dtype=bool)
    '''
    if not L:
        L=1
        while gt > ((1<<L) - 1):
            L += 1
    return _np.array(map(lambda l: bool(gt&(1<<(L-l-1))),range(L)))


def integerify(b):
    '''Transform a binary sequence on the HC into an integer.

    Parameters:
       - b: bool vector representing a genotype

    Returns:
       - gt: integer representing the same genotype

    **Examples**:

    .. sourcecode:: ipython

       In [1]: integerify([False, True, True])
       Out[1]: 3
    '''
    L = len(b)
    a = [(1<<(L-l-1)) for l in range(L)]
    return _np.dot(b,a)
%}
/*****************************************************************************/

/*****************************************************************************/
/* HAPLOID_LOWD                                                              */
/*****************************************************************************/
%define DOCSTRING_HAPLOID_LOWD
"Class for low-dimensional population genetics (short genomes ~20 loci).

The class offers a number of functions, but an example will explain the basic idea::

    #####################################
    #   EXAMPLE SCRIPT                  #
    #####################################
    import numpy as np
    import matplotlib.pyplot as plt
    import FFPopSim as h

    c = h.haploid_lowd(5)               # 5 loci

    # initialize with 300 individuals with genotype 00000,
    # and 700 with genotype 00010
    c.set_genotypes([0, 2], [300, 700])

    # set an additive fitness landscape with these coefficients
    c.set_fitness_additive([0.02,0.03,0.04,0.02, -0.03])
    # Note: we are in the -/+ basis, so
    #        F[10000] - F[00000] = 2 * 0.02
    # Hence the coefficients are half of the effect of mutation on fitness

    c.evolve(100)                       # evolve for 100 generations
    c.plot_diversity_histogram()
    plt.show()
    #####################################
"
%enddef
%feature("autodoc", DOCSTRING_HAPLOID_LOWD) haploid_lowd;
%extend haploid_lowd {
/* constructor */
%define DOCSTRING_HAPLOID_LOWD_INIT
"Construct a low-dimensional population with certain parameters.

Parameters:
    - L : number of loci (at least 1)
    - rng_seed : seed for the random number generator
"
%enddef
%feature("autodoc", DOCSTRING_HAPLOID_LOWD_INIT) haploid_lowd;
%exception haploid_lowd {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Construction impossible. Please check input args.");
                SWIG_fail;
        }
}

/* string representations */
%feature("autodoc", "x.__str__() <==> str(x)") __str__;
%feature("autodoc", "x.__repr__() <==> repr(x)") __repr__;
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"haploid_lowd: L = %d, N = %f", (int)$self->L(), $self->N());
        return &buffer[0];
}
const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"<haploid_lowd(%d, %g)>", (int)$self->L(), $self->N());
        return &buffer[0];
}

/* copy */
%pythoncode
%{
def copy(self, rng_seed=0):
    '''Copy population into new instance.

    Parameters:
       - rng_seed: random number to initialize the new population
    '''
    pop = haploid_lowd(self.L, rng_seed=rng_seed)

    # Mutation and recombination
    if self.recombination_model not in ['FREE_RECOMBINATION']:
        pop.set_recombination_rates(self.get_recombination_rates(), self.recombination_model)
    tmp = self.get_mutation_rates()
    if _np.isscalar(tmp) or tmp.ndim < 2:
        pop.set_mutation_rates(tmp)
    else:
        pop.set_mutation_rates(*tmp)
    pop.circular = self.circular
    pop.outcrossing_rate = self.outcrossing_rate

    # Fitness
    pop.set_fitness_function(range(1<<self.L), self.get_fitnesses())

    # Population parameters
    pop.carrying_capacity = self.carrying_capacity
    pop.set_genotypes(range(1<<self.L), self.get_genotype_frequencies() * self.N)

    # Evolution
    pop.generation = self.generation

    return pop
%}

/* ignore hypercubes */
%ignore fitness;
%ignore population;

/* read/write attributes */
%feature("autodoc", "is the genome circular?") circular;
%feature("autodoc", "current carrying capacity of the environment") carrying_capacity;
%feature("autodoc", "outcrossing rate (probability of sexual reproduction per generation)") outcrossing_rate;

/* read only attributes */
%ignore get_number_of_loci;
%feature("autodoc", "Number of loci (read-only)") L;
%feature("autodoc", "Number of loci (read-only)") number_of_loci;
const int L;
const int number_of_loci;

%ignore get_population_size;
%feature("autodoc", "Population size (read-only)") N;
%feature("autodoc", "Population size (read-only)") population_size;
const int N;
const int population_size;

%ignore get_generation;
%ignore set_generation;
%feature("autodoc", "Current generation (read-only)") generation;
int generation;

/* recombination model */
%rename (_get_recombination_model) get_recombination_model;
%rename (_set_recombination_model) set_recombination_model;
%feature("autodoc",
"Model of recombination to use

Available values:

   - FFPopSim.FREE_RECOMBINATION: free shuffling between parents
   - FFPopSim.CROSSOVERS: block recombination with crossover probability
   - FFPopSim.SINGLE_CROSSOVER: block recombination with crossover probability
") get_recombination_model;
%exception set_recombination_model {
  $action
  if (result == HG_BADARG) {
     PyErr_SetString(PyExc_ValueError,"Recombination model nor recognized.");
     SWIG_fail;
  } else if (result == HG_MEMERR) {
     PyErr_SetString(PyExc_ValueError,"Unable to allocate/release memory for the recombination patterns.");
     SWIG_fail;
  }
}
%pythoncode
%{
recombination_model = property(_get_recombination_model, _set_recombination_model)
%}

/* status function */
%pythoncode
%{
def status(self):
    '''Print a status list of the population parameters'''
    parameters = (('number of loci', 'L'),
                  ('circular', 'circular'),
                  ('population size', 'N'),
                  ('carrying capacity', 'carrying_capacity'),
                  ('generation', 'generation'),
                  ('outcrossing rate', 'outcrossing_rate'),
                  ('recombination model', 'recombination_model'),
                 )
    lenmax = max(map(lambda x: len(x[0]), parameters))

    for (strin, name) in parameters:
        par = getattr(self, name)
        # Recombination model needs a conversion
        # (a very frequently used one, to be honest)
        if strin == 'recombination model':
            if par == 0:
                par = 'FREE_RECOMBINATION'
            elif par == 1:
                par = 'SINGLE_CROSSOVER'
            else:
                par = 'CROSSOVERS'
        print(('{:<'+str(lenmax + 2)+'s}').format(strin)+'\t'+str(par))
%}

/* initialize frequencies */
%feature("autodoc",
"Initialize the population in linkage equilibrium with specified allele frequencies.

Parameters:
   - frequencies: an array of length L with all allele frequencies
   - N: set the population size and, if still unset, the carrying
     capacity to this value

.. note:: the population size is only used for resampling and has therefore
          no effect on the speed of the simulation.
") set_allele_frequencies;
%pythonprepend set_allele_frequencies {
if len(args) and (len(args[0]) != self.L):
    raise ValueError('The input array of allele frequencies has the wrong length.')
}
%exception set_allele_frequencies {
  $action
  if (result) {
     PyErr_SetString(PyExc_RuntimeError,"Error in the C++ function.");
     SWIG_fail;
  }
}
%pythonappend set_allele_frequencies {
return None
}

/* initialize genotypes */
%ignore set_genotypes;
%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* indices), (int len2, double* vals)};
int _set_genotypes(int len1, double* indices, int len2, double* vals) {
        vector<index_value_pair_t> gt;
        index_value_pair_t temp;
        for(size_t i = 0; i != (size_t)len1; i++) {
                temp.index = (int)indices[i];
                temp.val = vals[i];
                gt.push_back(temp);
        }
        return $self->set_genotypes(gt);
}
%clear (int len1, double* indices);
%clear (int len2, double* vals);
%pythoncode
%{
def set_genotypes(self, genotypes, counts):
    '''Initialize population with fixed counts for specific genotypes.

    Parameters:
       - genotypes: list of genotypes to set. Genotypes are specified as integers,
                    from 00...0 that is 0, up to 11...1 that is 2^L-1.
       - counts: list of counts for those genotypes

    .. note:: the population size and, if unset, the carrying capacity will be set as the sum of the counts.
    .. note:: you can use Python binary notation for the indices, e.g. 0b0110 is 6.
    '''
    genotypes = _np.asarray(genotypes, float)
    counts = _np.asarray(counts, float)
    if len(genotypes) != len(counts):
        raise ValueError('Indices and counts must have the same length')
    if self._set_genotypes(genotypes, counts):
        raise RuntimeError('Error in the C++ function.')
%}

/* initialize wildtype */
%feature("autodoc",
"Initialize population of N individuals with the - allele at all loci (wildtype)

Parameters:
   - N: the number of individuals

.. note:: the carrying capacity is set to the same value if still unset.
") set_wildtype;
%exception set_wildtype {
  $action
  if (result) {
     PyErr_SetString(PyExc_RuntimeError,"Error in the C++ function.");
     SWIG_fail;
  }
}
%pythonappend set_wildtype {
return None
}

/* recombination rates */
%feature("autodoc",
"Get the recombination rate between the specified locus and the following one.
") get_recombination_rate;
%pythonprepend get_recombination_rate {
if len(args) and (args[0] >= self.L - 1):
    raise ValueError("Expecting a locus from 0 to L - 2.")
}

%pythoncode
%{
def get_recombination_rates(self):
    '''Get recombination rates.

Returns:
    - the rates between neighboring loci, a list of float of length L-1

.. note:: if the recombination model if FREE_RECOMBINATION, an error is raised.
    '''
    if self.L < 2:
        raise ValueError('There is no recombination with less than 2 loci.')

    rm = self.recombination_model
    if rm == FREE_RECOMBINATION:
        raise ValueError(('The current recombination model is free recombination,'+
                          'hence recombination rates are not defined.'+
                          ' Could you possibly mean outcrossing rate?'))
    elif rm in [SINGLE_CROSSOVER, CROSSOVERS]:
        return _np.array([self.get_recombination_rate(i) for i in range(self.L - 1)])
    else:
        raise RuntimeError('Recombination model not found')
%}

%rename (_set_recombination_rates) set_recombination_rates;
%pythoncode
%{
def set_recombination_rates(self, rates, model=None):
    '''Set the recombination rate(s).

Parameters:
    - rates: if a double, the recombination rate at between any two loci; if an array,
      the locus-specific recombination rates
    - model: the recombination model to use (CROSSOVERS or, for linear
      genomes, SINGLE_CROSSOVER)

.. note:: if locus-specific rates are specified, the array must have length
          (L-1) for linear chromosomes and length L for circular ones. The
          i-th element is the crossover rate between the i-th site and the
          (i+1)-th site.

.. note:: if the recombination model is not specified, the current model will be kept or,
          if the current model is FREE_RECOMBINATION, then CROSSOVERS will be set.
    '''

    # Default recombination model
    if model is None:
        if self.recombination_model != FREE_RECOMBINATION:
            model = self.recombination_model
        else:
            model = CROSSOVERS

    # Check whether the model makes sense
    if model == FREE_RECOMBINATION:
        raise ValueError("Cannot assign rates to free recombination!")
    if model not in (CROSSOVERS, SINGLE_CROSSOVER):
        raise ValueError("Model not recognized.")
    if (self.circular and (model == SINGLE_CROSSOVER)):
        raise ValueError("Single crossover not available for circular genomes.")

    # Check whether the chromosome is circular
    if self.circular:
        len_rates = self.L
    else:
        len_rates = self.L - 1

    # Check whether the input argument is a list or a scalar
    if _np.isscalar(rates):
        self._set_recombination_rates([rates] * len_rates, model)

    elif len(rates) != len_rates:
        raise ValueError("Expecting an array of length "+str(len_rates)+".")
    else:
        self._set_recombination_rates(rates, model)

%}

/* mutation rate(s) */
%rename (_get_mutation_rate) get_mutation_rate;
%pythoncode
%{
def get_mutation_rates(self, locus=None, direction=None):
    '''Get one or several mutation rates.

Parameters:
    - locus: get only the mutation rate(s) of this locus
    - direction: get only the forward or backward mutation rate(s). This argument
                 is a Boolean, 0/False for forward rates, 1/True for backward rates.

Returns:
    - the mutation rate(s) requested

**Note**: if the mutation rates for all loci and/or directions are the same,
this function will try to be smart and give you the answer you are looking for.
In case of doubt, you will get a matrix (L x 2) with the full mutation rate
landscape.
    '''

    if locus is not None:
        if not _np.isscalar(locus):
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
            mrs = _np.array([self._get_mutation_rate(l, direction) for l in range(self.L)])
            if len(_np.unique(mrs)) == 1:
                return mrs[0]
            else:
                return mrs
        else:
            mrs = _np.array([[self._get_mutation_rate(l, d) for l in range(self.L)] for d in [0,1]])
            if len(_np.unique(mrs)) == 1:
                return mrs[0,0]
            else:
                return mrs
%}

%ignore set_mutation_rates;
int _set_mutation_rates(double *IN_ARRAY2, int DIM1, int DIM2) {
        double ** mrs = new double*[2];
        for(size_t i = 0; i < 2; i++)
                mrs[i] = &(IN_ARRAY2[DIM2 * i]);
        int result = $self->set_mutation_rates(mrs);
        delete[] mrs;
        return result;
}
%pythoncode
%{
def set_mutation_rates(self, rates, rates_back=None):
    '''Set the mutation rate(s).

Parameters:
    - rates:if a double, the mutation rate at any locus in both directions
      or, if rates_back is not None, only in the forward direction

      if a vector, the mutation rate is specified for each locus, the same
      in both directions or, if rates_back is not None, only in the
      forward direction

    - rates_back: mutation rate in the backward direction (global or
      locus-specific)
    '''

    L = self.L
    if _np.isscalar(rates):
        if rates_back is None:
            ratesm = _np.repeat(rates, L * 2).reshape(2,L)
        else:
            ratesm = _np.vstack([_np.repeat(rates, L), _np.repeat(rates_back, L)])
    elif (_np.rank(rates) != 1) or ((rates_back is not None) and (_np.rank(rates_back) != 1)):
        raise ValueError('Please input one/two numbers or arrays.')
    else:
        if rates_back is None:
            ratesm = _np.vstack([rates, rates])
        else:
            ratesm = _np.vstack([rates, rates_back])

    if self._set_mutation_rates(ratesm):
        raise RuntimeError('Error in the C++ function.')
%}

/* evolve */
%feature("autodoc",
"Evolve for some generations

Parameters:
    - gen: number of generations to evolve the population, defaults to one
") evolve;
%exception evolve {
  $action
  if (result) {
     PyErr_SetString(PyExc_RuntimeError,"Error in the C++ function.");
     SWIG_fail;
  }
}
%pythonappend evolve {
return None
}

%feature("autodoc",
"Evolve for some generations deterministically (skips the resampling)

Parameters:
    - gen: number of generations to evolve the population
") evolve_deterministic;
%exception evolve_deterministic {
  $action
  if (result) {
     PyErr_SetString(PyExc_RuntimeError,"Error in the C++ function.");
     SWIG_fail;
  }
}
%pythonappend evolve_deterministic {
return None
}

%feature("autodoc",
"Evolve for some generations without recombination

Parameters:
    - gen: number of generations to evolve the population
") evolve_norec;
%exception evolve_norec {
  $action
  if (result) {
     PyErr_SetString(PyExc_RuntimeError,"Error in the C++ function.");
     SWIG_fail;
  }
}
%pythonappend evolve_norec {
return None
}

/* get genotype frequencies */
%feature("autodoc",
"Get the frequency of a genotype

Parameters:
    - genotype: genotype, whose the frequency is to be returned

Returns:
    - the frequency of the genotype
") get_genotype_frequency;
%pythonprepend get_genotype_frequency {
if len(args) and (args[0] >= (1<<self.L)):
    raise ValueError("Expecting an individual from 0 to 2^L - 1.")
}

%pythoncode
%{
def get_genotype_frequencies(self):
    '''Get the frequency of each genotype.'''
    return _np.array([self.get_genotype_frequency(l) for l in range(1<<self.L)])
%}

/* get allele frequencies */
%feature("autodoc",
"Get the frequency of the + allele

Parameters:
    - locus: locus, at which the frequency of the + allele is to be computed

Returns:
    - the frequency of the + allele, :math:`\\nu_i := \\frac{1 + \\left<s_i\\right>}{2}`, where :math:`s_i \in \{-1, 1\}`.
") get_allele_frequency;
%pythonprepend get_allele_frequency {
if len(args) and (args[0] >= (self.L)):
    raise ValueError("Expecting a locus from 0 to L - 1.")
}

%pythoncode
%{
def get_allele_frequencies(self):
    '''Get the frequencies of all + alleles'''
    return _np.array([self.get_allele_frequency(l) for l in range(self.L)])
%}

%feature("autodoc",
"Get the frequency of genotypes with the + allele at both loci.

Parameters:
    - locus1: first locus
    - locus2: second locus

Returns:
    - the joint frequency of the + alleles
") get_pair_frequency;
%pythonprepend get_pair_frequency {
if (len(args) >= 2) and ((args[0] >= (self.L)) or (args[1] >= (self.L))):
    raise ValueError("Expecting loci from 0 to L - 1.")
}

%feature("autodoc",
"Get chi of an allele in the -/+ basis

Parameters:
    - locus: locus whose chi is to be computed

Returns:
    - the chi of that allele, :math:`\\chi_i := \\left<s_i\\right>`, where :math:`s_i \in \{-1, 1\}`.
") get_chi;
%pythonprepend get_chi {
if len(args) and (args[0] >= (self.L)):
    raise ValueError("Expecting a locus from 0 to L - 1.")
}

%feature("autodoc",
"Get :math:`\\chi_{ij}`

Parameters:
    - locus1: first locus
    - locus2: second locus

Returns:
    - the linkage disequilibiurm between them, i.e. :math:`\\chi_{ij} := \\left<s_i s_j\\right> - \\chi_i \\cdot \\chi_j`.
") get_chi2;
%pythonprepend get_chi2 {
if (len(args) >= 2) and ((args[0] >= (self.L)) or (args[1] >= (self.L))):
    raise ValueError("Expecting loci from 0 to L - 1.")
}

%feature("autodoc",
"Get linkage disequilibrium

Parameters:
    - locus1: first locus
    - locus2: second locus

Returns:
    - the linkage disequilibiurm between them, i.e. :math:`D_{ij} := 1 / 4 \\left[\\left<s_i s_j\\right> - \\chi_i \\cdot \\chi_j\\right]`.
") get_LD;
%pythonprepend get_LD {
if (len(args) >= 2) and ((args[0] >= (self.L)) or (args[1] >= (self.L))):
    raise ValueError("Expecting loci from 0 to L - 1.")
}

%feature("autodoc",
"Get moment of two alleles in the -/+ basis

Parameters:
    - locus1: first locus
    - locus2: second locus

Returns:
    - the second moment, i.e. :math:`\\left<s_i s_j\\right>`, where :math:`s_i, s_j \in \{-1, 1\}`.
") get_moment;
%pythonprepend get_moment {
if (len(args) >= 2) and ((args[0] >= (self.L)) or (args[1] >= (self.L))):
    raise ValueError("Expecting loci from 0 to L - 1.")
}

/* random sampling */
%pythoncode
%{
def random_genomes(self, n_sample):
    '''Get random genomes according sampled from the population.

    Parameters:
        - n_sample: number of random genomes to sample

    Returns:
        - integers corresponding to random genomes in the population.
    '''
    counts = _np.random.multinomial(n_sample, self.get_genotype_frequencies())
    ind = counts.nonzero()[0]
    counts = counts[ind]
    sample = _np.concatenate([_np.repeat(ind[i], counts[i]) for i in range(len(ind))])
    _np.random.shuffle(sample)
    return sample
%}

/* get fitness */
%feature("autodoc",
"Get fitness values of a genotype

Parameters:
    - genotype: genotype whose fitness is to be calculated. This can either be an integer or in binary format, e.g. 5 = 0b101

Returns:
    - the fitness of that genotype.
") get_fitness;
%pythonprepend get_fitness {
if len(args) and (args[0] >= (1<<self.L)):
    raise ValueError("Expecting an individual between 0 and 2^L - 1.")
}

%feature("autodoc", "Get the fitness of all possible genotypes.") get_fitnesses;
%pythonprepend get_fitnesses {
args = tuple([1<<self.L] + list(args))
}
void get_fitnesses(int DIM1, double* ARGOUT_ARRAY1) {
        for(size_t i=0; i < (size_t)DIM1; i++)
                ARGOUT_ARRAY1[i] = $self->get_fitness(i);
}

/* get fitness coefficients (Fourier space) */
%feature("autodoc",
"Get fitness coefficient of a combination (bitset) of loci

Parameters:
    - loci_bitset: Bitset of loci interested by the coefficient (see below). This can either be an integer or in binary format, e.g. 5 = 0b101

.. note:: Examples for loci_bitset:
   - 0: fitness baseline for the population
   - (1 << X): additive coefficient for locus X
   - (1 << X) + (1 << Y): epistatic (2-locus) coefficient between loci X and Y

.. note:: Remember that FFPopSim uses 1/-1 based hypercubes.

Returns:
    - the fitness coefficient of that combination of loci.
") get_fitness_coefficient;
%pythonprepend get_fitness_coefficient {
if len(args) and (args[0] >= (1<<self.L)):
    raise ValueError("Expecting an individual between 0 and 2^L - 1.")
}

%feature("autodoc",
"Get all fitness coefficients.

The order of the coefficient is by bitset of interested loci:
- the population baseline is at position 0
- the additive term of locus X is at position (1 << X), i.e. 2^X
- the 2-locus epistatic term between loci X and Y is at (1 << X) + (1 << Y), i.e. 2^X + 2^Y
and so on.

For instance, the following indices contain:
0 aka 0b0: baseline for the population
1 aka 0b1: additive coefficient for the first locus
2 aka 0b10: additive coefficient for the second locus
3 aka 0b11: epistatic coefficient bewteen locus one and two
4 aka 0b100: additive coefficient for the third locus
5 aka 0b101: epistatic coefficient between locus one and three
6 aka 0b110: epistatic coefficient between locus two and three
7 aka 0b111: epistatic coefficient among loci one, two, and three (3-locus term)
and so on.
") get_fitness_coefficients;
%pythonprepend get_fitness_coefficients {
args = tuple([1<<self.L] + list(args))
}
void get_fitness_coefficients(int DIM1, double* ARGOUT_ARRAY1) {
        for(size_t i=0; i < (size_t)DIM1; i++)
                ARGOUT_ARRAY1[i] = $self->get_fitness_coefficient(i);
}


/* divergence/diversity/fitness distributions and plot (full Python implementations) */
%pythoncode
%{
def get_fitness_histogram(self, n_sample=1000, **kwargs):
    '''Get the histogram of the fitness of a sample from the population.

    Parameters:
        - n_sample: number of individual to sample at random from the population. defaults to 1000

    Returns:
       - h: numpy.histogram of fitness in the population
    '''

    # Random sample
    gt = self.random_genomes(n_sample)

    # Calculate fitness
    fit = _np.array([self.get_fitness(gt[i]) for i in range(n_sample)])

    return _np.histogram(fit, **kwargs)


def plot_fitness_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot the histogram of the fitness of a sample from the population.

    Parameters:
        - axis: use an already existing axis for the plot
        - n_sample: number of individual to sample at random from the population. Defaults to 1000.
        - kwargs: further optional keyword arguments to numpy.histograms
    '''

    import matplotlib.pyplot as plt

    # Random sample
    gt = self.random_genomes(n_sample)

    # Calculate fitness
    fit = _np.array([self.get_fitness(gt[i]) for i in range(n_sample)])

    # Plot
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Fitness histogram')
        axis.set_xlabel('Fitness')
    axis.hist(fit, **kwargs)


def get_divergence_statistics(self, n_sample=1000):
    '''Get the mean and variance of the divergence of a population sample -- same as mean and variance of allele frequencies.

    Parameters:
        - n_sample: number of individuals to sample at random from the population. defaults to 1000.

    Returns:
        - stat: structure with mean and variance of divergence in the population
    '''

    L = self.L

    # Random sample
    gt = self.random_genomes(n_sample)

    # Calculate divegence
    div = _np.array([binarify(gt[i], L).sum() for i in range(n_sample)], int)

    return stat(div.mean(), div.var())


def get_divergence_histogram(self, bins=10, n_sample=1000, **kwargs):
    '''Get the histogram of the divergence of a population sample.

    Parameters:
        - bins: number of bins or list of bin edges (passed verbatim to numpy.histogram)
        - n_sample: number of individual to sample at random from the population, defaults to 1000.
        - kwargs: further optional keyword arguments to numpy.histograms

    Returns:
       - h: numpy.histogram of divergence in the population

    *Note*: to get a normalized histogram, use the *density* keyword.
    '''

    # Random sample
    gt = self.random_genomes(n_sample)

    # Calculate divergence
    div = _np.array([binarify(gt[i], self.L).sum() for i in range(n_sample)], int)

    return _np.histogram(div, bins=bins, **kwargs)


def plot_divergence_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot the histogram of the divergence of a population sample.

    Parameters:
        - axis: use an already existing axis for the plot
        - n_sample: number of individual to sample at random from the population, defaults to 1000.
        - kwargs: further optional keyword arguments to numpy.histograms
    '''
    import matplotlib.pyplot as plt
    L = self.L

    # Random sample
    gt = self.random_genomes(n_sample)

    # Calculate divegence
    div = _np.array([binarify(gt[i], L).sum() for i in range(n_sample)], int)

    # Plot
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Divergence histogram')
        axis.set_xlabel('Divergence')

    if 'bins' not in kwargs:
        kwargs['bins'] = _np.arange(10) * max(1, (div.max() + 1 - div.min()) / 10) + div.min()
    axis.hist(div, **kwargs)


def get_diversity_statistics(self, n_sample=1000):
    '''Get the mean and variance of the diversity of a population sample

    Parameters:
        - n_sample: number of individual to sample at random from the population, defaults to 1000.

    Returns:
        - stat: structure with mean and variance of diversity in the population
    '''

    # Random sample
    gt1 = self.random_genomes(n_sample)
    gt2 = self.random_genomes(n_sample)

    # Calculate diversity
    div = _np.array([binarify(gt1[i] ^ gt2[i], self.L).sum() for i in range(n_sample)], int)

    return stat(div.mean(), div.var())


def get_diversity_histogram(self, bins=10, n_sample=1000, **kwargs):
    '''Get the histogram of the diversity in a sample from the population.

    Parameters:
        - bins: number of bins or list of bin edges (passed verbatim to numpy.histogram)
        - n_sample: number of individual to sample at random from the population, defaults to 1000.
        - kwargs: further optional keyword arguments to numpy.histograms

    Returns:
       - h: numpy.histogram of diversity in the population

    *Note*: to get a normalized histogram, use the *density* keyword.
    '''

    # Random sample
    gt1 = self.random_genomes(n_sample)
    gt2 = self.random_genomes(n_sample)

    # Calculate diversity
    div = _np.array([binarify(gt1[i] ^ gt2[i], self.L).sum() for i in range(n_sample)], int)

    # Calculate histogram
    return _np.histogram(div, bins=bins, **kwargs)


def plot_diversity_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot the histogram of the diversity of a population sample.

    Parameters:
        - axis: use an already existing axis for the plot
        - n_sample: number of individual to sample at random from the population, defaults to 1000.
        - kwargs: further optional keyword arguments to numpy.histograms
    '''
    import matplotlib.pyplot as plt

    # Random sample
    gt1 = self.random_genomes(n_sample)
    gt2 = self.random_genomes(n_sample)

    # Calculate diversity
    div = _np.array([binarify(gt1[i] ^ gt2[i], self.L).sum() for i in range(n_sample)], int)

    # Plot
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Diversity histogram')
        axis.set_xlabel('Diversity')

    if 'bins' not in kwargs:
        kwargs['bins'] = _np.arange(10) * max(1, (div.max() + 1 - div.min()) / 10) + div.min()
    axis.hist(div, **kwargs)
%}

/* set fitness landscape */
%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* indices), (int len2, double* vals)};
int _set_fitness_func(int len1, double* indices, int len2, double* vals) {
        vector<index_value_pair_t> iv;
        index_value_pair_t temp;
        for(size_t i = 0; i != (size_t)len1; i++) {
                temp.index = (int)indices[i];
                temp.val = vals[i];
                iv.push_back(temp);
        }
        return ($self->fitness).init_list(iv);
}

int _set_fitness_coeff(int len1, double* indices, int len2, double* vals) {
        vector<index_value_pair_t> iv;
        index_value_pair_t temp;
        for(size_t i = 0; i != (size_t)len1; i++) {
                temp.index = (int)indices[i];
                temp.val = vals[i];
                iv.push_back(temp);
        }
        return ($self->fitness).init_coeff_list(iv);
}
%clear (int len1, double* indices);
%clear (int len2, double* vals);
%pythoncode
%{
def set_fitness_function(self, genotypes, values):
    '''Set the fitness landscape for individual genotypes.

    Parameters:
       - genotypes: genotype to which the fitness values will be assigned. Genotypes are specified as integers,
                    from 00...0 that is 0, up to 11...1 that is 2^L-1.
       - values: fitness values to assign

    .. note:: you can use Python binary notation for the genotypes, e.g. 0b0110 is 6.
    '''
    genotypes = _np.asarray(genotypes, float)
    values = _np.asarray(values, float)
    if len(genotypes) != len(values):
        raise ValueError('Indices and values must have the same length')
    if self._set_fitness_func(genotypes, values):
        raise RuntimeError('Error in the C++ function.')


def set_fitness_coefficients(self, coefficients, values):
    '''Set the fitness landscape in Fourier space for individual Fourier coefficients.

    Parameters:
       - coefficients: Fourier coefficients to which the values will be assigned. They are specified
                       as integers, from 00...0 that is 0, up to 11...1 that is 2^L-1.
       - values: values to assign

    .. note:: you can use Python binary notation for the coefficients, e.g. 0b0110 is 6.
    '''
    coefficients = _np.asarray(coefficients, float)
    values = _np.asarray(values, float)
    if len(coefficients) != len(values):
        raise ValueError('Indices and values must have the same length')
    if self._set_fitness_coeff(coefficients, values):
        raise RuntimeError('Error in the C++ function.')

%}

/* set additive fitness component */
%feature("autodoc",
"Set an additive fitness landscape. Coefficients obey +/- convention.

Parameters:
    - coefficients: array/list of additive fitness coefficients. It must have length L.
") set_fitness_additive;
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

/* entropy */
%feature("autodoc",
"Get the genotype entropy of the population

.. note:: the genotype entropy is defined as :math:`-\\sum_{i=0}^{2^L} p_i \\log p_i`.
") genotype_entropy;
%feature("autodoc",
"get the allele entropy of the population

.. note:: the allele entropy is defined as :math:`-\\sum_{i=0}^{L} \\left[\\nu_i\log \\nu_i + (1-\\nu_i)\log(1-\\nu_i)\\right]`.
") allele_entropy;

/* ignore tests */
%ignore test_recombinant_distribution();
%ignore test_recombination(double *rec_rates);
%ignore mutation_drift_equilibrium(double** mutrates);
} /* extend haploid_lowd */

%{
const int haploid_lowd_L_get(haploid_lowd *h) {
  return (const int) h->get_number_of_loci();
}
const int haploid_lowd_number_of_loci_get(haploid_lowd *h) {
  return (const int) h->get_number_of_loci();
}

const int haploid_lowd_N_get(haploid_lowd *h) {
  return (const int) h->get_population_size();
}
const int haploid_lowd_population_size_get(haploid_lowd *h) {
  return (const int) h->get_population_size();
}

int haploid_lowd_generation_get(haploid_lowd *h) {
  return (const int) h->get_generation();
}
void haploid_lowd_generation_set(haploid_lowd *h, int g) {
  h->set_generation(g);
}
%}
/*****************************************************************************/
