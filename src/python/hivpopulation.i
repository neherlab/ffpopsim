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

/*****************************************************************************/
/* HIVGENE                                                                   */
/*****************************************************************************/
%feature("autodoc", "Structure for an HIV gene.") hivgene;
%extend hivgene {
%feature("autodoc", "x.__str__() <==> str(x)") __str__;
%feature("autodoc", "x.__repr__() <==> repr(x)") __repr__;
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"hivgene: start: %d, end: %d", $self->start, $self->end);
        return &buffer[0];
}
const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"hivgene(%d, %d)", $self->start, $self->end);
        return &buffer[0];
}

%feature("autodoc", "Initial position of the gene") start;
%feature("autodoc", "Final position of the gene") end;
}

/*****************************************************************************/
/* HIVPOPULATION                                                             */
/*****************************************************************************/
%define DOCSTRING_HIVPOPULATION
"Class for HIV population genetics (genome size = 10000).

This class is the main object for simulating the evolution of HIV.
The class offers a number of functions, but an example will explain the basic
idea::

   #####################################
   #   EXAMPLE SCRIPT                  #
   #####################################
   import numpy as np
   import matplotlib.pyplot as plt
   import FFPopSim as h
   
   c = h.hivpopulation(2000)        # Create a population of 2000 individuals
   c.evolve(100)                    # Evolve (neutrally) for 100 generations
   c.plot_divergence_histogram()
   plt.show()
   #####################################

**This class is a subclass of haploid_high and offers most of its methods.**
In addition to the haploid_highd class, this class offers functions for reading
fitness and drug resistance landscapes from a text file, and to save genomes as
plain text or in compressed NumPy format.

Moreover, there are two phenotypic traits, replication and resistance. Their
relative importance for viral fitness is set by the ``treatment`` attribute::

   f[trait] = trait[0] + treatment * trait[1]

By default, ``treatment`` is set to zero, to simulate non-treated patients.

The gene structure of HIV is not modelled explicitely, except for a stub of
1000 sites between position 7000 and 8000 to roughly model the _env_ gene.
"
%enddef
%feature("autodoc", DOCSTRING_HIVPOPULATION) hivpopulation;

%extend hivpopulation {
%feature("autodoc",
"Construct a HIV population with certain parameters.

Parameters:

   - N     number of viral particles
   - rng_seed	seed for the random number generator. If this is 0, time(NULL)+getpid() is used.
   - mutation_rate	mutation rate in events / generation / site
   - coinfection_rate	probability of coinfection of the same cell by two viral particles in events / generation
   - crossover_rate	probability of template switching during coinfection in events / site

.. note:: the genome length is 10000 (see HIVGENOME).
") hivpopulation;
/* constructor */
%exception hivpopulation {
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
        sprintf(buffer,"hivpopulation: N = %d", $self->N());
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"<hivpopulation(%d)>", $self->N());
        return &buffer[0];
}

/* copy */
%pythoncode{
def copy(self, rng_seed=0):
    '''Copy population into new instance.
    
    Parameters:
       - rng_seed: random number to initialize the new population
    '''
    pop = hivpopulation(self.N,
                        rng_seed=rng_seed,
                        mutation_rate=self.mutation_rate,
                        coinfection_rate=self.outcrossing_rate,
                        crossover_rate=self.crossover_rate)

    # Fitness
    for i in xrange(self.number_of_traits):
        pop.set_trait_additive(self.get_trait_additive(i), i)
        for coeff in self.get_trait_epistasis(i):
            pop.add_trait_coefficient(coeff[0], coeff[1], i)

    # Population parameters
    pop.carrying_capacity = self.carrying_capacity
    pop.set_genotypes(self.get_genotypes(), self.get_clone_sizes())    

    # Evolution
    pop._set_generation(self.generation)
    
    return pop
}

/* we have two traits anyway */
%ignore add_fitness_coefficient;
%ignore clear_fitness;

/* treatment */
%ignore set_treatment;
%ignore get_treatment;
%feature("autodoc",
"Treatment weight (between 0 and 1)

.. note:: this variable controls how important is either of the two phenotypic
          traits, replication and resistance. Their contribution to fitness is
          always linear (in this implementation).
") treatment;
double treatment;

/* read selection/resistance coefficients */
%feature("autodoc",
"Read replication coefficients from a text file

Parameters:
   - filename: string with the name of the file to read the coefficients from
") read_replication_coefficients;

%feature("autodoc",
"Read resistance coefficients from a text file

Parameters:
   - filename: string with the name of the file to read the coefficients from
") read_resistance_coefficients;

/* write genotypes */
%feature("autodoc",
"Store random genotypes into a plain text file.

Parameters:
   - filename: string with the name of the file to store the genotype into
   - sample_size: how many random genotypes to store
   - gt_label: common fasta label for the genotypes (e.g. 'HIV-sim')
   - start: if only a portion of the genome is to be stored, start from this position
   - length: store a chunk from ``start`` to this length
") write_genotypes;

%pythoncode {
def write_genotypes_compressed(self, filename, sample_size, gt_label='', start=0, length=0):
    '''Store random genotypes into a compressed file.

    Parameters:
       - filename: string with the name of the file to store the genotype into
       - sample_size: how many random genotypes to store
       - gt_label: common fasta label for the genotypes (e.g. "HIV-sim")
       - start: if only a portion of the genome is to be stored, start from this position
       - length: store a chunk from ``start`` to this length

    The genotypes can be read using numpy.load.
    '''

    import numpy as np 
    L = self.number_of_loci
    if length <= 0:
        length = L - start
    d = {}
    for i in xrange(sample_size):
        rcl = self.random_clone()
        d['>'+str(i)+'_GT-'+gt_label+'_'+str(rcl)] = self.get_genotype(rcl)[start:start+length]
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
    '''Set HIV trait landscape according to some general parameters.

    Parameters:
       - lethal_fraction: fraction of lethal sites
       - deleterious_fraction: fraction of deleterious sites
       - adaptive_fraction: fraction of beneficial sites
       - effect_size_lethal: effect of lethal changes
       - effect_size_deleterious: average effect of deleterious changes
       - effect_size_adaptive: average effect of beneficial changes

       - env_fraction: fraction of beneficial sites in env
       - effect_size_env: average effect of beneficial changes in env
       - number_epitopes: number of (epistatic) epitopes
       - epitope_strength: average height of an epitope escape mutation
       - number_valleys: number of (epistatic) valleys
       - valley_strength: average depth of a valley

    .. note:: the effects of deleterious and beneficial sites are exponentially
              distributed, i.e. most of them will still be almost neutral.
    
    .. note:: fractions refer to first and second positions only. For instance,
              by default, 80% of first and second positions outside env are
              deleterious.

    .. note:: the third positions are always neutral (synonymous).
    '''

    import numpy as np
    
    # Clear trait
    self.clear_trait(traitnumber)

    # Handy
    L = self.L
    aL = np.arange(L)

    # Decide what mutation is of what kind
    # Note: the rest, between
    #
    # lethal_fraction + deleterious_fraction and (1 - adaptive_fraction),
    #
    # is neutral, i.e. EXACTLY 0. Fair assumption.
    onetwo_vector = (aL % 3) < 2
    random_numbers = np.random.random(L)
    adaptive_mutations = (random_numbers > (1 - adaptive_fraction)) & onetwo_vector
    lethal_mutations = (random_numbers < lethal_fraction) & onetwo_vector
    deleterious_mutations = ((random_numbers > lethal_fraction) & \
                             (random_numbers < (lethal_fraction + deleterious_fraction)) & \
                             (random_numbers < (1 - adaptive_fraction)) & \
                             onetwo_vector)
    
    # Decide how strong mutations are
    single_locus_effects=np.zeros(L)
    single_locus_effects[np.where(deleterious_mutations)] = -np.random.exponential(effect_size_deleterious, deleterious_mutations.sum())
    single_locus_effects[np.where(adaptive_mutations)] = np.random.exponential(effect_size_adaptive, adaptive_mutations.sum())
    single_locus_effects[np.where(lethal_mutations)] = -effect_size_lethal
    
    # Mutations in env are treated separately
    env_position = (aL >= self.env.start) & (aL < self.env.end)
    env_mutations = (random_numbers > (1 - env_fraction)) & onetwo_vector & env_position
    single_locus_effects[np.where(env_mutations)] = np.random.exponential(effect_size_env, env_mutations.sum())
        
    # Call the C++ routines
    self.set_trait_additive(single_locus_effects, traitnumber)

    # Epistasis
    multi_locus_coefficients=[]
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

    # Set fitness valleys
    for vi in xrange(number_valleys):
        pos = np.random.random_integers(L/3-100)
        d = int(np.random.exponential(10) + 1)
        valley_str = np.random.exponential(valley_strength)
        if number_valleys:
            print 'valley:', pos*3, valley_str
        (f1,f2,f12)=add_valley(valley_str)
        single_locus_effects[pos*3+1]+=f1
        single_locus_effects[(pos+d)*3+1]+=f2
        multi_locus_coefficients.append([[pos*3+1, (pos+d)*3+1], f12])
    
    # Set epitopes (bumps, i.e. f_DM < d_WT << f_SM)
    for ei in xrange(number_epitopes):
        pos = np.random.random_integers(L/3-10)
        epi_strength = np.random.exponential(epitope_strength)
        if number_epitopes:
                print 'epitope', pos*3, epi_strength
        epi, f1,f2,f12=add_epitope(epi_strength)
        single_locus_effects[(pos+epi[0])*3+1]+=f1
        single_locus_effects[(pos+epi[1])*3+1]+=f2
        multi_locus_coefficients.append([[(pos+epi[0])*3+1, (pos+epi[1])*3+1], f12])

    for mlc in multi_locus_coefficients:
        self.add_trait_coefficient(mlc[1], np.asarray(mlc[0], int), traitnumber)
    self._update_traits()
    self._update_fitness()
}

/* helper functions for replication and resistance */
/* There is a reason why they are not properties, namely because you will never be able
to set them by slicing, e.g. pop.additive_replication[4:6] = 3. In order to implement
this functionality we would need a whole subclass of ndarray with its own set/get
methods, and nobody is really keen on doing this. */
%pythoncode{
def get_replication_additive(self):
    '''The additive part of the replication lansdscape.

    Returns:
       - coefficients: array of additive replication coefficients

    .. warning:: the -/+ basis is used throughout the library.
                 If you are used to the 0/1 basis, keep in mind that
                 the interaction series-expansion is different.
    '''
    return self.get_trait_additive(0)


def set_replication_additive(self, coefficients):
    '''Set the additive replication coefficients

    Parameters:
       - coefficients: array of additive replication coefficients

    .. warning:: the -/+ basis is used throughout the library.
                 If you are used to the 0/1 basis, keep in mind that
                 the interaction series-expansion is different.
    '''

    self.set_trait_additive(coefficients, 0)


def get_resistance_additive(self):
    '''The additive part of the resistance lansdscape.

    Returns:
       - coefficients: array of additive drug resistance coefficients

    .. warning:: the -/+ basis is used throughout the library.
                 If you are used to the 0/1 basis, keep in mind that
                 the interaction series-expansion is different.
    '''
    return self.get_trait_additive(1)


def set_resistance_additive(self, coefficients):
    '''Set the additive drug resistance coefficients

    Parameters:
       - coefficients: array of additive drug resistance coefficients

    .. warning:: the -/+ basis is used throughout the library.
                 If you are used to the 0/1 basis, keep in mind that
                 the interaction series-expansion is different.
    '''

    self.set_trait_additive(coefficients, 1)


}

/* Generate random landscapes */
%pythoncode{
def set_replication_landscape(self,
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
    '''Set the phenotypic landscape for the replication capacity of HIV.
    
    Parameters:
       - lethal_fraction: fraction of lethal sites
       - deleterious_fraction: fraction of deleterious sites
       - adaptive_fraction: fraction of beneficial sites
       - effect_size_lethal: effect of lethal changes
       - effect_size_deleterious: average effect of deleterious changes
       - effect_size_adaptive: average effect of beneficial changes

       - env_fraction: fraction of beneficial sites in env
       - effect_size_env: average effect of beneficial changes in env
       - number_epitopes: number of (epistatic) epitopes
       - epitope_strength: average height of an epitope escape mutation
       - number_valleys: number of (epistatic) valleys
       - valley_strength: average depth of a valley

    .. note:: the effects of deleterious and beneficial sites are exponentially
              distributed, i.e. most of them will still be almost neutral.
    
    .. note:: fractions refer to first and second positions only. For instance,
              by default, 80% of first and second positions outside env are
              deleterious.

    .. note:: the third positions are always neutral (synonymous).
    '''

    self.set_trait_landscape(traitnumber=0,
                        lethal_fraction=lethal_fraction,
                        deleterious_fraction=deleterious_fraction,
                        adaptive_fraction=adaptive_fraction,
                        effect_size_lethal=effect_size_lethal,
                        effect_size_deleterious=effect_size_deleterious,
                        effect_size_adaptive=effect_size_adaptive,
                        env_fraction=env_fraction,
                        effect_size_env=effect_size_env,
                        number_epitopes=number_epitopes,
                        epitope_strength=epitope_strength,
                        number_valleys=number_valleys,
                        valley_strength=valley_strength)


def set_resistance_landscape(self,
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
    '''Set the phenotypic landscape for the drug resistance of HIV.
    
    Parameters:
       - lethal_fraction: fraction of lethal sites
       - deleterious_fraction: fraction of deleterious sites
       - adaptive_fraction: fraction of beneficial sites
       - effect_size_lethal: effect of lethal changes
       - effect_size_deleterious: average effect of deleterious changes
       - effect_size_adaptive: average effect of beneficial changes

       - env_fraction: fraction of beneficial sites in env
       - effect_size_env: average effect of beneficial changes in env
       - number_epitopes: number of (epistatic) epitopes
       - epitope_strength: average height of an epitope escape mutation
       - number_valleys: number of (epistatic) valleys
       - valley_strength: average depth of a valley

    .. note:: the effects of deleterious and beneficial sites are exponentially
              distributed, i.e. most of them will still be almost neutral.
    
    .. note:: fractions refer to first and second positions only. For instance,
              by default, 80% of first and second positions outside env are
              deleterious.

    .. note:: the third positions are always neutral (synonymous).
    '''

    self.set_trait_landscape(traitnumber=1,
                        lethal_fraction=lethal_fraction,
                        deleterious_fraction=deleterious_fraction,
                        adaptive_fraction=adaptive_fraction,
                        effect_size_lethal=effect_size_lethal,
                        effect_size_deleterious=effect_size_deleterious,
                        effect_size_adaptive=effect_size_adaptive,
                        env_fraction=env_fraction,
                        effect_size_env=effect_size_env,
                        number_epitopes=number_epitopes,
                        epitope_strength=epitope_strength,
                        number_valleys=number_valleys,
                        valley_strength=valley_strength)


}
} /* extend hivpopulation */

%{
double hivpopulation_treatment_get(hivpopulation *h) {
  return (double) h->get_treatment();
}
void hivpopulation_treatment_set(hivpopulation *h, double t) {
  h->set_treatment(t);
}
%} /* attributes of hivpopulation */
/*****************************************************************************/
