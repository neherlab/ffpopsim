/**** HIVGENE ****/
%define DOCSTRING_HIVGENE
"Structure for an HIV gene."
%enddef
%feature("autodoc", DOCSTRING_HIVGENE) hivgene;

%extend hivgene {
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"start: %d, end: %d", $self->start, $self->end);
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"hivgene(%d, %d)", $self->start, $self->end);
        return &buffer[0];
}
}


/**** HIVPOPULATION ****/
%define DOCSTRING_HIVPOPULATION
"Class for HIV population genetics (genome size = 10000).

This class is the main object for simulating the evolution of HIV.
The class offers a number of functions, but an example will explain the basic
idea:

#####################################
#   EXAMPLE SCRIPT                  #
#####################################
import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h

c = h.hivpopulation(2000)
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

In addition to the haploid_highd class, this class offers functions for reading
fitness and drug resistance landscapes from a text file, and to save genomes as
plain text or in compressed numerical Python format.
"
%enddef
%feature("autodoc", DOCSTRING_HIVPOPULATION) hivpopulation;

%extend hivpopulation {

%define DOCSTRING_HIVPOPULATION_INIT
"Construct a HIV population with certain parameters.

Parameters:
- N     number of viral particles
- rng_seed	seed for the random number generator. If this is 0, time(NULL)+getpid() is used.
- mutation_rate	mutation rate in events / generation / site
- coinfection_rate	probability of coinfection of the same cell by two viral particles in events / generation
- crossover_rate	probability of template switching during coinfection in events / site

Note: the genome length is 10000 (see HIVGENOME).
"
%enddef
%feature("autodoc", DOCSTRING_HIVPOPULATION_INIT) hivpopulation;

/* we have two traits anyway */
%ignore add_fitness_coefficient;
%ignore clear_fitness;

/* constructor */
%exception hivpopulation {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Construction impossible. Please check input args.");
                SWIG_fail;
        }
}

/* treatment */
%rename (_set_treatment) set_treatment;
%rename (_get_treatment) get_treatment;
%pythoncode {
treatment = property(_get_treatment, _set_treatment)
}

/* read selection/resistance coefficients */
%typemap(in) istream &model (std::ifstream temp) {
        if (!PyString_Check($input)) {
                PyErr_SetString(PyExc_ValueError, "Expecting a string");
                return NULL;
        }
        temp.open(PyString_AsString($input));
        $1 = &temp;
}

/* write genotypes */
%typemap(in) ostream &out_genotypes (std::ofstream temp) {
        if (!PyString_Check($input)) {
                PyErr_SetString(PyExc_ValueError, "Expecting a string");
                return NULL;
        }
        temp.open(PyString_AsString($input));
        $1 = &temp;
}

%pythoncode {
def write_genotypes_compressed(self, filename, sample_size, gt_label='', start=0, length=0):
        '''Write genotypes into a compressed archive.'''
        import numpy as np 
        L = self.number_of_loci
        if length <= 0:
                length = L - start
        d = {}
        for i in xrange(sample_size):
                rcl = self.random_clone()
                d['>'+str(i)+'_GT-'+gt_label+'_'+str(rcl)] = self.get_genotype(rcl,L)[start:start+length]
        np.savez_compressed(filename, **d)    
}


/* set trait landscape */
%pythoncode {
def set_trait_landscape(self,
                        traitnumber=0,
                        lethal_fraction=0.05,
                        deleterious_fraction=0.8,
                        adaptive_fraction=0.01,
                        effect_size_lethal=0.8,
                        effect_size_deleterious=0.1,
                        effect_size_adaptive=0.01,
                        env_fraction=0.1,
                        effect_size_env=0.01,
                        number_epitopes=0,
                        epitope_strength=0.05,
                        number_valleys=0,
                        valley_strength=0.1,
                        ):
    '''Set HIV trait landscape according to some general parameters.'''

    import numpy as np
    from scipy import stats

    # Nested functions
    def add_epitope(strength=0.2):
        '''Note: we are in the +-1 basis.'''
        loci = random.sample(range(9),2)
        loci.sort()
        depression = - 0.05
        f1 = depression*0.25
        f2 = depression*0.25
        f12 = depression*0.25 - strength*0.5
        return loci, f1,f2,f12
     
    def add_valley(depth=0.1, height=0.01):
        '''Note: we are in the +-1 basis.'''
        f1 = height*0.25
        f2 = height*0.25
        f12 = height*0.25 + depth*0.5
        return (f1,f2,f12)

    
    L = self.L
    aL = np.arange(L)

    # Initialize fitness coefficients as zero (neutral model)
    single_locus_effects=np.zeros(L)
    multi_locus_coefficients=[]
            
    # Set single locus fitness coefficients
    first_codon_position = np.arange(0,L,3)
    second_codon_position = np.arange(1,L,3)
    # Note: the third positions are always neutral (synonymous)
    
    # Decide what mutation is of what kind
    onetwo_vector = (aL % 3) < 2
    random_numbers = np.random.random(L)
    adaptive_mutations = (random_numbers > (1 - adaptive_fraction)) & onetwo_vector
    lethal_mutations = (random_numbers < lethal_fraction) & onetwo_vector
    deleterious_mutations = (random_numbers > lethal_fraction) & (random_numbers < (lethal_fraction + deleterious_fraction)) & onetwo_vector
    
    # Decide how strong mutations are
    adaptive_dis = stats.expon(scale=effect_size_adaptive)
    deleterious_dis = stats.expon(scale=effect_size_deleterious)
    single_locus_effects[np.where(adaptive_mutations)] += adaptive_dis.rvs(adaptive_mutations.sum())
    single_locus_effects[np.where(deleterious_mutations)] -= deleterious_dis.rvs(deleterious_mutations.sum())
    single_locus_effects[np.where(lethal_mutations)] -= effect_size_lethal
    
    # Mutations in env are treated separately
    env_position = (aL >= self.env.start) * (aL < self.env.end)
    env_mutations = (random_numbers>(1 - env_fraction)) & onetwo_vector & env_position
    env_dis = stats.expon(scale=effect_size_env)
    single_locus_effects[np.where(env_mutations)] += env_dis.rvs(env_mutations.sum())
        
    # Note: the rest, between lethal_fraction + deleterious_fraction and (1 -
    # adaptive_fraction), is neutral, i.e. EXACTLY 0. Fair assumption.
    
    # Set fitness valleys
    for vi in xrange(number_valleys):
        pos = np.random.random_integers(L/3-100)
        d = int(stats.expon(scale=10).rvs() +1)
        valley_dis=stats.expon(scale=valley_strength)
        valley_str = valley_dis.rvs()
        if number_valleys:
            print 'valley:', pos*3, valley_str
        (f1,f2,f12)=add_valley(valley_str)
        single_locus_effects[pos*3+1]+=f1
        single_locus_effects[(pos+d)*3+1]+=f2
        multi_locus_coefficients.append([[pos*3+1, (pos+d)*3+1], f12])
    
    # Set epitopes (bumps, i.e. f_DM < d_WT << f_SM)
    for ei in xrange(number_epitopes):
        pos = np.random.random_integers(L/3-10)
        epi_dis=stats.expon(scale=epitope_strength)
        epi_strength = epi_dis.rvs()
        if number_epitopes:
                print 'epitope', pos*3, epi_strength
        epi, f1,f2,f12=add_epitope(epi_strength)
        single_locus_effects[(pos+epi[0])*3+1]+=f1
        single_locus_effects[(pos+epi[1])*3+1]+=f2
        multi_locus_coefficients.append([[(pos+epi[0])*3+1, (pos+epi[1])*3+1], f12])
                

    # Call the C++ routines
    self.clear_trait(traitnumber)
    self.set_additive_trait(single_locus_effects, traitnumber)
    for mlc in multi_locus_coefficients:
        self.add_trait_coefficient(mlc[1], np.asarray(mlc[0], int), traitnumber)
    self.calc_stat()
}

/* helper functions for replication and resistance */
%pythoncode{
def get_additive_replication(self):
        '''Get the additive part of the replication lansdscape.'''
        return self.get_additive_trait(0)


def get_additive_resistance(self):
        '''Get the additive part of the resistance lansdscape.'''
        return self.get_additive_trait(1)


def set_additive_replication(self, single_locus_effects):
        '''Set the additive part of the replication lansdscape.'''
        self.set_additive_trait(single_locus_effects, 0)


def set_additive_resistance(self, single_locus_effects):
        '''Set the additive part of the resistance lansdscape.'''
        self.set_additive_trait(single_locus_effects, 1)


def set_replication_landscape(self, **kwargs):
        '''Set the phenotypic landscape for the replication capacity of HIV.
        
        Parameters:
        -  traitnumber=0
        -  lethal_fraction=0.05
        -  deleterious_fraction=0.8
        -  adaptive_fraction=0.01
        -  effect_size_lethal=0.8
        -  effect_size_deleterious=0.1
        -  effect_size_adaptive=0.01
        -  env_fraction=0.1
        -  effect_size_env=0.01
        -  number_epitopes=0
        -  epitope_strength=0.05
        -  number_valleys=0
        -  valley_strength=0.1
        '''
        kwargs['traitnumber']=0
        self.set_trait_landscape(**kwargs)


def set_resistance_landscape(self, **kwargs):
        '''Set the phenotypic landscape for the drug resistance of HIV.
        
        Parameters:
        -  traitnumber=0
        -  lethal_fraction=0.05
        -  deleterious_fraction=0.8
        -  adaptive_fraction=0.01
        -  effect_size_lethal=0.8
        -  effect_size_deleterious=0.1
        -  effect_size_adaptive=0.01
        -  env_fraction=0.1
        -  effect_size_env=0.01
        -  number_epitopes=0
        -  epitope_strength=0.05
        -  number_valleys=0
        -  valley_strength=0.1 
        '''

        kwargs['traitnumber']=0
        self.set_trait_landscape(**kwargs)
}

} /* extend hivpopulation */
