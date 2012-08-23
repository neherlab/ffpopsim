/* renames and ignores */
%ignore coeff_t;
%ignore coeff_single_locus_t;
%ignore hypercube_highd;

/* general typemaps */
/* convert a Python bool array into a boost::dynamic_bitset */
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

/**** CLONE_T ****/
%feature("autodoc", "Clone of isogenetic individuals") clone_t;

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

/* read/write attributes */
%feature("autodoc", "Fitness of the clone") fitness;
%feature("autodoc", "Number of individuals in this clone") clone_size;

/* traits */
%feature("autodoc", "Number of traits (read-only)") _get_number_of_traits;
%rename (_trait) trait;
int _get_number_of_traits() {
        return ($self->trait).size();
}
%pythoncode {
number_of_traits = property(_get_number_of_traits)
}

void _get_trait(int DIM1, double* ARGOUT_ARRAY1) {
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = ($self->trait)[i];
}
%pythoncode {
@property
def trait(self):
    '''Traits vector of the clone'''
    return self._get_trait(self.number_of_traits)
}

/* genotype */
%ignore genotype;
int _get_genotype_length() {return ($self->genotype).size();}
void _get_genotype(int DIM1, short* ARGOUT_ARRAY1) {
        for(size_t i=0; i < ($self->genotype).size(); i++) ARGOUT_ARRAY1[i] = ($self->genotype)[i];
}

void _set_genotype(boost::dynamic_bitset<> genotype_in) {$self->genotype = genotype_in;}

%pythoncode {
@property
def genotype(self):
    '''Genotype of the clone'''
    import numpy as np
    return np.array(self._get_genotype(self._get_genotype_length()), bool)


@genotype.setter
def genotype(self, genotype):
    self._set_genotype(genotype)
}
} /* extend clone_t */

/**** HAPLOID_HIGHD ****/
%define DOCSTRING_HAPLOID_HIGHD
"Class for high-dimensional population genetics (genomes larger than ~20 loci).

This class is the main object for simulating the evolution of populations with
many loci (more than ~20). The class offers a number of functions, but an
example will explain the basic idea::

   ######################################
   #  EXAMPLE SCRIPT FOR HAPLOID_HIGHD  #
   ######################################
   import numpy as np
   import matplotlib.pyplot as plt
   import FFPopSim as h
   c = h.haploid_highd(300)       # 300 loci
   pop.set_wildtype(1000)         # start with 1000 wildtype individuals
   pop.mutation_rate = 1e-7       # mutation rate per site per generation
   pop.outcrossing_rate = 1e-1    # probability of sexual reproduction per gen
   pop.crossover_rate = 1e-2      # probability of crossover per site per gen
   pop.evolve(100)                # evolve for 100 generations
   c.plot_divergence_histogram()
   plt.show()
   ######################################

Populations can have a number of phenotypic traits that concur to the fitness
of each individual. The function that calculates fitness from the phenotype
identifies fitness with the first trait only by default. The user is, however,
free to subclass haploid_highd in C++ (as it is done in hivpopulation) and
implement their own phenotype -> fitness function.

In addition, the trait landscapes describe the genotype -> phenotype maps.
These can be set directly from Python (since the genotypic space has a finite
number of elements).

**Note**: fitness is not a phenotypic trait directly, but rather a function of *all*
phenotypic traits together. 
"
%enddef
%feature("autodoc", DOCSTRING_HAPLOID_HIGHD) haploid_highd;

%extend haploid_highd {

/* constructor */
%feature("autodoc",
"Construct a high-dimensional population with certain parameters.

Parameters:
- L     length of the genome(number of loci)
- rng_seed      seed for the random generator. If zero (default) pick a random number
- number_of_traits      number of phenotypic traits
") haploid_highd;

/* TODO: ignore hypercubes for now */
%ignore trait;

/* TODO: ignore pointers to the clones for now */
%ignore current_pop;
%ignore new_pop;

/* ignore stream stuff we never need in Python */
%ignore print_allele_frequencies;
%ignore get_genotype_string;
%ignore read_ms_sample;
%ignore read_ms_sample_sparse;

/* ignore weird functions using pointers */
%ignore get_pair_frequencies(vector < vector <int> > *loci);
%ignore random_clones(unsigned int n_o_individuals, vector <int> *sample);

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

/* read/write attributes */
%feature("autodoc", "is the genome circular?") circular;
%feature("autodoc", "current carrying capacity of the environment") carrying_capacity;
%feature("autodoc", "outcrossing rate (probability of sexual reproduction per generation)") outcrossing_rate;
%feature("autodoc", "crossover rate (probability of crossover per site per generation)") crossover_rate;
%feature("autodoc", "mutation rate (per site per generation)") mutation_rate;
%feature("autodoc",
"model of recombination to use

Available values:
   - FFPopSim.FREE_RECOMBINATION: free shuffling between parents
   - FFPopSim.CROSSOVERS: block recombination with crossover probability
") recombination_model;


/* read only parameters */
%ignore L;
%ignore N;
%rename (_get_number_of_loci) get_number_of_loci;
%rename (_get_population_size) get_population_size;
%rename (_get_generation) get_generation;
%rename (_get_number_of_clones) get_number_of_clones;
%rename (_get_number_of_traits) get_number_of_traits;
%rename (_get_max_fitness) get_max_fitness;
%rename (_get_participation_ratio) get_participation_ratio;
%pythoncode {
L = property(_get_number_of_loci)
N = property(_get_population_size)
number_of_loci = property(_get_number_of_loci)
population_size = property(_get_population_size)
generation = property(_get_generation)
number_of_clones = property(_get_number_of_clones)
number_of_traits = property(_get_number_of_traits)
max_fitness = property(_get_max_fitness)
participation_ratio = property(_get_participation_ratio)
}
%feature("autodoc", "Number of loci (read-only)") get_number_of_loci;
%feature("autodoc", "Population size (read-only)") get_population_size;
%feature("autodoc", "Number of clones (read-only)") get_number_of_clones;
%feature("autodoc", "Number of traits (read-only)") get_number_of_traits;
%feature("autodoc", "Current generation (read-only)") get_generation;
%feature("autodoc", "Maximal fitness in the population (read-only)") get_max_fitness;
%feature("autodoc", "Participation ratio (read-only)") get_participation_ratio;

/* initialize wildtype */
%feature("autodoc",
"Set up a population of wildtype individuals

Parameters:
   - N: the number of individuals and carrying capacity
") set_wildtype;

/* initalize frequencies */
%ignore set_allele_frequencies;
int _set_allele_frequencies(double *IN_ARRAY1, int DIM1, int n_o_genotypes) {return $self->set_allele_frequencies(IN_ARRAY1, n_o_genotypes);}
%pythoncode {
def set_allele_frequencies(self, frequencies, N):
    '''Initialize the population according to the given allele frequencies.

    Parameters:
       - frequencies: an array of length L with all allele frequencies
       - N: the number of individuals and carrying capacity
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
       - indices: list of genotypes to set (e.g. 0 -> 00...0, L-1 -> 11...1)
       - counts: list of counts for those genotypes

    .. note:: the population size and the carrying capacity are set as the sum of the counts.
    .. note:: you can use Python binary notation for the indices, e.g. 0b0110 = 6.

    **Example**: if you want to initialize 200 individuals with genotype 001 and 300 individuals
                 with genotype 110, you can use ``set_genotypes([[0,0,1], [1,1,0]], [200, 300])``
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

/* bottleneck */
%feature("autodoc",
"Make the population undergo a bottleneck

Parameters:
   - size_of_bottleneck: the number of individuals at the bottleneck
") bottleneck;

/* statistics */
%feature("autodoc", "Calculate trait and fitness statistics for the population") calc_stat;

%feature("autodoc",
"Get the mean and variance of the divergence in the population.

Parameters:
   - n_sample: number of individuals to sample at random from the population

Returns:
   - stat: structure with mean and variance of divergence in the population
") get_divergence_statistics;

%feature("autodoc", 
"Get the mean and variance of the diversity in the population.

Parameters:
   - n_sample: number of individuals to sample at random from the population

Returns:
   - stat: structure with mean and variance of diversity in the population
") get_diversity_statistics;

%feature("autodoc",
"Get the mean and variance of a trait in the population.

Parameters:
   - t: number of the trait whose statistics are to be computed

Returns:
   - stat: structure with mean and variance of the trait in the population
") get_trait_statistics;

%feature("autodoc",
"Get the mean and variance of the fitness in the population.

Returns:
   - stat: structure with mean and variance of the fitness in the population
") get_fitness_statistics;

%feature("autodoc",
"Get the covariance of two traits in the population.

Parameters:
   - t1: first trait
   - t2: second trait

Returns:
   - cov: the covariance of the two traits
") get_trait_covariance;



/* get allele frequencies */
void _get_allele_frequencies(double* ARGOUT_ARRAY1, int DIM1) {
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = $self->get_allele_frequency(i);
}
%pythoncode {
def get_allele_frequencies(self):
    '''Get all allele frequencies'''
    return self._get_allele_frequencies(self.L)
}

%feature("autodoc",
"Get the frequency of the + allele at the selected locus

Parameters:
   - locus: locus whose frequency of the + allele is to be computed

Returns:
   - frequency: allele frequency in the population
") get_allele_frequency;

%feature("autodoc",
"Get chi of an allele in the -/+ basis

Parameters:
    - locus: locus whose chi is to be computed

Returns:
    - the chi of that allele, :math:`\\chi_i := \\left<s_i\\right>`, where :math:`s_i \in \{-1, 1\}`.
") get_chi;

%feature("autodoc",
"Get the joint frequency of two + alleles

Parameters:
    - locus1: first locus
    - locus2: second locus

Returns:
    - the joint frequency of the + alleles
") get_pair_frequency;


/* get genotypes */
void _get_genotype(unsigned int n, short* ARGOUT_ARRAY1, int DIM1) {
        boost::dynamic_bitset<> newgt = (*($self->current_pop))[n].genotype;
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = newgt[i];
}
%pythoncode {
def get_genotype(self, n):
    '''Get a genotype from the population

    Parameters:
       - n: index of the clone whose genotype is to be returned

    Returns:
       - gt: boolean array of the genotype
    '''

    return self._get_genotype(n, self.number_of_loci)


def get_genotypes(self, ind=None):
    '''Get genotypes of the population.

    Parameters:
       - ind: if a scalar, a single genotype corresponding to clone ind is returned;
         otherwise, several genotypes are returned (default: all)
    '''

    import numpy as np
    L = self.number_of_loci
    if np.isscalar(ind):
        return np.array(self._get_genotype(ind, L), bool)

    if ind is None:
        ind = xrange(self.number_of_clones)
    genotypes = np.zeros((len(ind), L), bool)
    for i, indi in enumerate(ind):
        genotypes[i] = self._get_genotype(indi, L)
    return genotypes
}

/* add genotypes */
%feature("autodoc",
"Add new individuals to the population with certain genotypes

Parameters:
   - gt: genotype to add to the population
   - n: number of new individuals carrying that genotype

.. note:: gt is an array/list that must be convertible into a bool.
") add_genotypes;

/* set trait/fitness coefficients */
%exception clear_fitness {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Fitness depends only on traits, not on the genome directly.");
                SWIG_fail;
        }
}

/* get single locus effects */
void _get_additive_trait(double* ARGOUT_ARRAY1, int DIM1, int t) {
        /* check trait number */
        if(t >= $self->get_number_of_traits())
                throw HP_BADARG;

        /* Initialize to zero */
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = 0;

        /* Add any coefficient you found */
        hypercube_highd *trait = &(($self->trait)[t]);
        coeff_single_locus_t * coeff;
        for(size_t i=0; i < trait->coefficients_single_locus.size(); i++) {
                coeff = &(trait->coefficients_single_locus[i]);
                ARGOUT_ARRAY1[coeff->locus] += coeff->value;
        }
}
%pythoncode{
def get_additive_trait(self, t=0):
    '''Get an array with the additive part of a trait for all loci.

    Parameters:
       - t: number of the trait

    Returns:
       - coefficients: array of additive coefficients for the selected trait
    '''
    return self._get_additive_trait(self.L, t)
}

/* update functions we should not need from Python */
%rename (_update_traits) update_traits;
%rename (_update_fitness) update_fitness;

/* set single locus effects */
%feature("autodoc",
"Set the additive part of a trait

Parameters:
   - coefficients: array of coefficients for the trait (of length L)
   - t: number of the trait to set
") set_additive_trait;
void set_additive_trait(int DIM1, double* IN_ARRAY1, int t=0) {
        /* check trait number */
        if(t >= $self->get_number_of_traits())
                throw HP_BADARG;
        /* check length of vector */
        if(DIM1 != $self->L())
                throw HP_BADARG;

        /* reset trait landscape */
        $self->trait[t].reset_additive();
        
        /* set the new coefficients */
        vector <int> loci(1,0);
        for(size_t i = 0; i < DIM1; i++) {
                if(abs(IN_ARRAY1[i]) > HP_NOTHING) {
                        loci[0] = i;
                        $self->add_trait_coefficient(IN_ARRAY1[i], loci, t);
                }
        }

        /* update the population */
        $self->update_traits();
        $self->update_fitness();
}

%feature("autodoc", "Shortcut for set_additive_trait when there is only one trait") set_additive_fitness;
void set_additive_fitness(int DIM1, double *IN_ARRAY1) {
        /* check whether we really have only one trait */
        if($self->get_number_of_traits() > 1)
                throw HP_BADARG;
        /* check length of vector */
        if(DIM1 != $self->L())
                throw HP_BADARG;

        /* reset trait landscape */
        $self->trait[0].reset_additive();
        
        /* set the new coefficients */
        vector <int> loci(1,0);
        for(size_t i = 0; i < DIM1; i++) {
                if(abs(IN_ARRAY1[i]) > HP_NOTHING) {
                        loci[0] = i;
                        $self->add_trait_coefficient(IN_ARRAY1[i], loci, 0);
                }
        }

        /* update the population */
        $self->update_traits();
        $self->update_fitness();
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
%feature("autodoc",
"Add a coefficient to the trait landscape.

Parameters:
   - value: value of the coefficient
   - loci: array/list of loci indexed by the coefficient.
   - t: number of the trait to be changed

**Example**: to set a second-order epistatic term :math:`t_{ij} = 0.1`, use ``add_trait_coefficient(0.1, [i, j])``.

.. warning:: the -/+ basis is used throughout the library. If you are used to the 0/1 basis, keep in mind that the interaction series-expansion is different.
") add_trait_coefficient;

%feature("autodoc", "Shortcut for add_trait_coefficient when there is only one trait") add_fitness_coefficient;

%feature("autodoc",
"Clear a trait landscape.

Parameters:
   - t: number of the trait to be cleared
") clear_trait;

%feature("autodoc", "Shortcut for clear_trait when there is only one trait") clear_fitness;
%feature("autodoc", "Clear all trait landscapes") clear_traits;

/* random epistasis */
%feature("autodoc",
"Set a random epistatic term in the genotype-phenotype map

Parameters:
   - epistasis_std: standard deviation of the random epistatic terms

.. note:: the epistatic terms will be Gaussian distributed around zero with the given standard deviation.
") set_random_trait_epistasis;

%feature("autodoc", "Shortcut for set_random_trait_epistasis when there is only one trait") set_random_epistasis;

/* fitness/traits of clones */
void _get_fitnesses(int DIM1, double* ARGOUT_ARRAY1) {
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = $self->get_fitness(i);
}
%pythoncode {
def get_fitnesses(self):
        '''Get the fitness of all clones.'''
        return self._get_fitnesses(self.number_of_clones)
}

%feature("autodoc",
"Get the fitness of an individual

Parameters:
   - n: index of the clone whose fitness is to be computed

Returns:
   - fitness: fitness value of that clone
") get_fitness;

%feature("autodoc",
"Get a trait of an individual

Parameters:
   - n: index of the clone whose trait is to be computed
   - t: trait to be computed

Returns:
   - trait: value of that trait for that clone
") get_trait;

/* get sizes of all clones */
void _get_clone_sizes(int DIM1, int* ARGOUT_ARRAY1) {
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = $self->get_clone_size(i);
}
%pythoncode {
def get_clone_sizes(self):
        '''Get the fitness of all clones.'''
        return self._get_clone_sizes(self.number_of_clones)
}

%feature("autodoc",
"Get the size of a clone

Parameters:
   - n: index of the clone

Returns:
   - size: size of the selected clone
") get_clone_size;

/* clone structure */
%feature("autodoc",
"Recompress the clone structure

During its evolution, identical clones might be generated by different routes at
different times. This function merges any such duplicates into unique clones with
the size equal to the sum of the sizes of the duplicates.
") unique_clones;

/* Hamming distance (full Python reimplementation) */
%ignore distance_Hamming;
%pythoncode {
def distance_Hamming(self, clone_gt1, clone_gt2, chunks=None, every=1):
    '''Calculate the Hamming distance between two genotypes

    Parameters:
       - clone_gt1: index of the clone corresponding to the first genotype
       - clone_gt2: index of the clone corresponding to the second genotype
       - chunks: list of pairs delimiting the genetic areas to include
       - every: do the comparison only on certain sites

    **Example**: to calculate the distance between the first two clones
    limited to third codon positions between locus 90 and 200, use:
    ``distance_Hamming(0, 1, chunks=[92, 200], every=3)``.
    '''
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
%pythoncode {
def random_genomes(self, n):
    '''Get random genomes from the population

    Parameters:
       - n: number of random genomes to compute

    Returns:
       - gts: (n x L) bool matrix with the n genotypes
    '''

    import numpy as np
    L = self.number_of_loci
    genotypes = np.zeros((n, L), bool)
    for i in xrange(genotypes.shape[0]):
        genotypes[i] = self._get_genotype(self.random_clone(), L)
    return genotypes
}

%feature("autodoc",
"Get random clones

Parameters:
   - n: number of random clones to return

Returns:
   - clones: clone indices
") random_clones;
void random_clones(int DIM1, unsigned int * ARGOUT_ARRAY1) {
        vector <int> sample = vector <int>(0);
        int err = $self->random_clones(DIM1, &sample);
        if(!err)
                for(size_t i=0; i < DIM1; i++)
                        ARGOUT_ARRAY1[i] = sample[i];
}

%feature("autodoc",
"Get a random clone

Returns:
   - clone: index of the random clone
") random_clone;

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

    Returns:
       - h: numpy.histogram of fitness in the population
    '''

    import numpy as np
    fit = [self.get_fitness(self.random_clone()) for i in xrange(n_sample)]
    h = np.histogram(fit, bins=bins, **kwargs)
    return h
    
    
def plot_fitness_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot a distribution of fitness in the population.

    Parameters:
       - axis: an axis to use. A new figure is created by default
       - n_sample: number of individuals to sample
       - kwargs: further optional keyword arguments to matplotlib.pyplot.hist

    Returns:
       - return value of axis.hist(...)
    '''

    import matplotlib.pyplot as plt
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
       - kwargs: further optional keyword arguments to numpy.histogram

    Returns:
       - h: numpy.histogram of divergence in the population
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
       - kwargs: further optional keyword arguments to matplotlib.pyplot.hist
    
    Returns:
       - return value of axis.hist(...)
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
       - kwargs: further optional keyword arguments to numpy.histogram

    Returns:
       - h: numpy.histogram of diversity in the population
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
       - kwargs: further optional keyword arguments to matplotlib.pyplot.hist
    
    Returns:
       - return value of axis.hist(...)
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
