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
import PopGenLib as h

c = h.hivpopulation(2000)
c.init_genotypes() 
c.evolve(10)
c.plot_divergence_histogram()
plt.show()
#####################################

An effective way to discover all available methods is to import PopGenLib from
an interactive shell (e.g. iPython), create a population as above, and use TAB
autocompletion:

In [1]: import PopGenLib as h
In [2]: c = h.haploid_highd(5000, 2000)
In [3]: c.      <--- TAB

In addition to the haploid_highd class, this class offers functions for reading
fitness and drug resistance landscapes from a text file, and to save genomes as
plain text or in compressed numerical Python format.
"
%enddef
%feature("autodoc", DOCSTRING_HIVPOPULATION) hivpopulation;

%extend hivpopulation {

/* we do not need this and it conflicts with base class constructor calling */
%ignore set_up;

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
@property
def treatment(self):
        return self._get_treatment()


@treatment.setter
def treatment(self, value):
        self._set_treatment(value)
}

/* read selection/resistance coefficients */
%ignore read_replication_coefficients;
%rename (read_replication_coefficients) _read_replication_coefficients;
int _read_replication_coefficients(char *model){
        ifstream modelstream(model);
        return $self->read_replication_coefficients(modelstream);
}

%ignore read_resistance_coefficients;
%rename (read_resistance_coefficients) _read_resistance_coefficients;
int _read_resistance_coefficients(char *model){
        ifstream modelstream(model);
        return $self->read_resistance_coefficients(modelstream);
}


/* write genotypes */
%ignore write_genotypes;
%rename (write_genotypes) _write_genotypes;
int _write_genotypes(char * filename, int sample_size, char * gt_label=NULL, int start=0, int length=0) {
        if(HIVPOP_VERBOSE >= 1) cerr<<"hivpython::write_genotypes(char * filename, int sample_size, char * gt_label, int start, int length)...";

        ofstream genotype_file(filename);        
        string gt_labelstring;
        if(gt_label == NULL)
                gt_labelstring = "";
        else
                gt_labelstring = string(gt_label);

        return $self->write_genotypes(genotype_file, sample_size, gt_labelstring, start, length);
}


%pythoncode {
def write_genotypes_compressed(self, filename, sample_size, gt_label='', start=0, length=0):
        '''Write genotypes into a compressed archive.'''
        import numpy as np 
        L = self.get_number_of_loci()
        if length <= 0:
                length = L - start
        d = {}
        for i in xrange(sample_size):
                rcl = self.random_clone()
                d['>'+str(i)+'_GT-'+gt_label+'_'+str(rcl)] = self._get_genotype(rcl,L)[start:start+length]
        np.savez_compressed(filename, **d)    
}


/* create trait (fitness) landscapes */
void clear_trait(unsigned int traitnumber) {
        if(traitnumber >= $self->get_number_of_traits())
                throw HIVPOP_BADARG;
        else
                ($self->trait)[traitnumber].reset();
}

/* single locus effects */
%apply (int DIM1, double* IN_ARRAY1) {(int L, double* single_locus_effects)};
void _set_additive_trait(int L, double* single_locus_effects, int traitnumber=0) {
        vector <int> loci;
        for(size_t i = 0; i < L; i++) {
                loci.push_back(i);
                $self->add_trait_coefficient(single_locus_effects[i], loci, traitnumber);
                loci.clear();
        }
}
%clear (int L, double* single_locus_effects);

/* glue code */
%pythoncode {
def _set_trait_landscape(self,
                        traitnumber=0,
                        lethal_fraction=0.05,
                        deleterious_fraction=0.8,
                        adaptive_fraction=0.01,
                        effect_size_lethal=0.01,
                        effect_size_deleterious=0.1,
                        effect_size_adaptive=0.8,
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

    
    L = self.L()
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
    self._set_additive_trait(single_locus_effects, traitnumber)
    for mlc in multi_locus_coefficients:
        self._add_trait_coefficient(mlc[1], np.asarray(mlc[0], int), traitnumber)


def set_replication_landscape(self, **kwargs):
        '''Set the phenotypic landscape for the replication capacity of HIV.
        
        Parameters:
        ---  traitnumber=0
        ---  lethal_fraction=0.05
        ---  deleterious_fraction=0.8
        ---  adaptive_fraction=0.01
        ---  effect_size_lethal=0.01
        ---  effect_size_deleterious=0.1
        ---  effect_size_adaptive=0.8
        ---  env_fraction=0.1
        ---  effect_size_env=0.01
        ---  number_epitopes=0
        ---  epitope_strength=0.05
        ---  number_valleys=0
        ---  valley_strength=0.1
        '''
        kwargs['traitnumber']=0
        self._set_trait_landscape(**kwargs)


def set_resistance_landscape(self, **kwargs):
        '''Set the phenotypic landscape for the drug resistance of HIV.
        
        Parameters:
        ---  traitnumber=0
        ---  lethal_fraction=0.05
        ---  deleterious_fraction=0.8
        ---  adaptive_fraction=0.01
        ---  effect_size_lethal=0.01
        ---  effect_size_deleterious=0.1
        ---  effect_size_adaptive=0.8
        ---  env_fraction=0.1
        ---  effect_size_env=0.01
        ---  number_epitopes=0
        ---  epitope_strength=0.05
        ---  number_valleys=0
        ---  valley_strength=0.1 
        '''

        kwargs['traitnumber']=0
        self._set_trait_landscape(**kwargs)
}

} /* extend hivpopulation */
