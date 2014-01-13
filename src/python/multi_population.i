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

%define DOCSTRING_MULTI_POPULATION
/* FIXME */   
"Class for high-dimensional population genetics of geographically divided populations with 
large genomes (genomes length more than ~20 loci).


This class is the main object for simulating the evolution of several geographically divided populations. 
The class offers a number of functions, but an example will explain the basic idea::

   ########################################
   #  EXAMPLE SCRIPT FOR MULTI_POPULATION #
   ########################################
   import numpy as np
   import matplotlib.pyplot as plt
   import FFPopSim as h
   c = h.multi_population(10, 300)  # 10 'locations' 300 loci
   pop.add_random_genotype(1000)    # start with 1000 individuals having the same random genotype
   pop.mutation_rate = 1e-4         # mutation rate per site per generation
   pop.outcrossing_rate = 1e-1      # probability of sexual reproduction per gen
   pop.crossover_rate = 1e-2        # probability of crossover per site per gen
   pop.evolve(100)                  # evolve for 100 generations
   c.plot_divergence_histogram()
   plt.show()
   ######################################

Populations can evolve in different gepgraphic 'locations' separately as well as with migration between them.
Each 'location' may be characterized by its trait->fitness function so that if the organisms migrate between them ,they may experience different evolutionary forces in different 'locations'.
The trait->fitness maps are defined at the beginning of the simuations but can be re-defined during the simulations. 

**Note**: fitness is not a phenotypic trait directly, but rather a function of *all*
phenotypic traits together. 
"
%enddef
%feature("autodoc", DOCSTRING_MULTI_POPULATION) multi_population;

%extend multi_population {
                
/* constructor */
/* FIXME */
%feature("autodoc",
"Construct a high-dimensional population with certain parameters.

   Parameters:
   - l: number of locations 
   - L: number of loci
   - rng_seed: seed for the random generator. If zero (default) pick a random number
   - number_of_traits: number of phenotypic traits, defaults to one
   - all_polymorphic: option to use an infinite-sites model tracking ancestral alleles
                      (only available with a single phenotypic trait and zero mutation rate)
") multi_population;
%exception multi_population {
        try {
                $action
        } catch (int err) {
            if (err == HP_BADARG)
                   PyErr_SetString(PyExc_ValueError,
                                "Construction impossible. Please check input args.");
            else if(err == HP_MEMERR)
                   PyErr_SetString(PyExc_MemoryError,
                                "Construction impossible. Out of memory.");

            SWIG_fail;
        }
}




/*Properties*/

%feature("autodoc",
"Set the maximal migrating part of the population.

Parameters:
   - critical_migration_rate: Maximal rate of the popualtion that can migrate to other locations.

") set_critical_migration_rate;
%exception set_critical_migration_rate {
    try {
        $action
    }catch (int err) {
        PyErr_SetString(PyExc_RuntimeError,"The value of the migration rate must be from 0 to 1");
        SWIG_fail;  
    }
}
%rename(_get_critical_migration_rate) get_critical_migration_rate;
%pythoncode{
@property
def critical_migration_rate(self):
   '''Critical migration rate (Read-only)'''
   return self._get_critical_migration_rate()
}

/*Migration rates*/
%ignore set_migration_rates (vector < vector < double > > new_migration_rates);

%feature("autodoc", 
"
Set matrix for the migration rates between locations.

Parameters: 
   - migration_rates: migration rates matrix. Rows of the matrix are source locations, columns - destinations.

Returns:
   - zero if successful
")set_migration_rates;
%pythonprepend set_migration_rates {
    m_rates = args[0]
    m_rates = _np.array(m_rates, float, copy=False, ndmin=2) 
    args = tuple ([m_rates.ravel()])
}
%exception set_migration_rates {
    try{
         $action
    }catch (int err) {
        PyErr_SetString(PyExc_RuntimeError,"Wrong migration rate values.\nDo they exceed the critical migration rate of the population?");          
        SWIG_fail;                                              
    }
}
%apply (int DIM1, double* IN_ARRAY1){(int len1, double* m_rates)};
int set_migration_rates(int len1, double* m_rates){
    vector < vector < double > > migration_rates;
    vector <double> temp_m_rates;
    for (int i = 0; i < $self->get_number_of_locations(); i ++){
        temp_m_rates.clear();
        for (int j = 0; j < $self->get_number_of_locations(); j ++){
            temp_m_rates.push_back(m_rates[i * $self->get_number_of_locations() + j]);
        }
        migration_rates.push_back(temp_m_rates);
    }
    return $self->set_migration_rates(migration_rates);
}
%clear (int len1, double* migration_rates);


%ignore get_number_of_locations;
%rename(_get_number_of_locations) get_number_of_locations;
%pythoncode{
@property
def number_of_locations(self):
    '''Number of locations (read-only)'''
    return self._get_number_of_locations()
}

%rename(_max_fitness) max_fitness;
%pythoncode{
@property
def max_fitness(self):
    '''Maximal fitnes of the population (read-only)'''
    return self._max_fitness()
}

%rename (_get_generation) get_generation;
%pythoncode{
@property
def generation(self):
    '''Current generation (read-only)'''
    return self._get_generation()
}

%rename(_number_of_migration_events) number_of_migration_events;
%pythoncode{
@property
def number_of_migrations(self):
    '''Number of the migration events (read-only)'''
    return self._number_of_migration_events
}

/*Carrying capacities*/
%feature("autodoc", 
"
Set location-specific carrying capacities of the population.

Parameters: 
   - carrying_capacities: list of carrying capacity values for each location

Returns:
   - zero if successful
")set_carrying_capacities;
%exception set_carrying_capacities {
        try{
            $action
        }  catch (int err) {
            if (err == HP_BADARG)
                   PyErr_SetString(PyExc_ValueError,
                                "Size mismatch between the input array and number of locations. \nYou must provide carrying capacity for each location once.");
            SWIG_fail;
        }
}

%ignore get_carrying_capacities;
void _get_carrying_capacities(int DIM1, double* ARGOUT_ARRAY1) {
        for(int i=0; i < (int)DIM1; i++)
                ARGOUT_ARRAY1[i] = ($self->get_carrying_capacity(i));
}
%pythoncode {
@property
def carrying_capacities(self):
    '''Location-specific carrying capacities of the population (read-only)'''
    return self._get_carrying_capacities(self.number_of_locations)
}

%ignore get_carrying_capacity;

/*Mutation rates*/
%feature("autodoc", 
"
Set location-specific mutation rates of the population.

Parameters: 
   - mutation_rates: list of mutation rate values for each location

Returns:
   - zero if successful
")set_mutaiton_rates;
%exception set_mutaiton_rates {
        try{
            $action
        }  catch (int err) {
            if (err == HP_BADARG)
                   PyErr_SetString(PyExc_ValueError,
                  "Size mismatch between the input array and number of locations. \nYou must provide mutaiton rate for each location once.");
            SWIG_fail;
        }
}

%ignore get_mutation_rates;
void _get_mutation_rates(int DIM1, double* ARGOUT_ARRAY1) {
        for(size_t i=0; i < (size_t)DIM1; i++)
                ARGOUT_ARRAY1[i] = ($self->get_mutation_rate(i));
}
%pythoncode {
@property
def mutation_rates(self):
    '''Location-specific mutation rates of the population (read-only)'''
    return self._get_mutation_rates(self.number_of_locations)
}



%ignore N;
void _N(int DIM1, int* ARGOUT_ARRAY1){
    for (int i = 0; i < DIM1; i++){
        ARGOUT_ARRAY1[i] = $self->N(i);
    } 
}
%pythoncode {
@property
def N(self):
   '''Population size (read-only)'''
   return self._N(self.number_of_locations)
}


%rename (_L) L;
%pythoncode{
@property
def L(self):
   '''Number of loci (read-only)'''
   return self._L()
}

/* trait weights */
%ignore set_trait_weights;
%pythonprepend _set_trait_weights {
if len(args) and (len(args[0]) != self.get_number_of_traits()):
    raise ValueError('The weights must be a sequence of length equal to the number of traits.')
}
void _set_trait_weights(double* IN_ARRAY1, int DIM1) {
        /* call the C++ method */
        $self->set_trait_weights(IN_ARRAY1);
        $self->update_fitness();
}
%feature("autodoc",
"weight of each trait on fitness

.. note:: Fitness is updated automatically when the weights are changed.
") _get_trait_weights;
%pythonprepend _get_trait_weights {
    args = tuple(list(args) + [self.get_number_of_traits()])
}
void _get_trait_weights(double* ARGOUT_ARRAY1, int DIM1) {
        /* check trait number */
        if(DIM1 != $self->get_number_of_traits())
                throw HP_BADARG; 

        /* set the output array */
        for(size_t t=0; t < (size_t)DIM1; t++)
                ARGOUT_ARRAY1[t] = $self->get_trait_weight(t);
}
%pythoncode {
trait_weights = property(_get_trait_weights, _set_trait_weights)
}

/* trait coefficients */
%feature("autodoc",
"Add a coefficient to the trait landscape.
 
Parameters:
   - value: value of the coefficient
   - loci: array/list of loci indexed by the coefficient.
   - t: number of the trait to be changed

**Example**: to set a second-order epistatic term :math:`t_{ij} = 0.1`, use ``add_trait_coefficient(0.1, [i, j])``.

.. warning:: the -/+ basis is used throughout the library. If you are used to the 0/1 basis, keep in mind that the interaction series-expansion is different.
") set_trait_coefficient;

/* genealogy */
%feature("autodoc",
"Track the genealogy of some loci.

   Parameters:
   - loci: sites whose genealogy is being stored

   Returns:
   - zero if successful, error code otherwise
") track_locus_genealogy;



/* implement multi_locus_genealogy as a read-only property */
%ignore genealogy;
%feature("autodoc",
"Genealogy of the tracked loci.

.. note:: This attribute is read-only.
") _get_genealogy;
multi_locus_genealogy _get_genealogy() {
        return $self->genealogy;
}
%pythoncode {
   genealogy = property(_get_genealogy)
}



/*evolution*/
%feature("autodoc",
"Evolve for some generations.

Parameters:
   - gen: number of generations, defaults to one
") evolve;
%exception evolve {
        try {
                $action
        } catch (int err) {
            if (err == HP_EXTINCTERR)
                   PyErr_SetString(PyExc_ValueError,
                                "Unfortunately, the popualtion went extinct. No further simulation possible. Aborting.");
            else 
                   PyErr_SetString(PyExc_ValueError,
                                "Unknown exception is thrown. Aborting.");

            SWIG_fail;
        }
}




/* initialize population */
%feature("autodoc",
"Initialize a population of N individuals with the same random genotype

Parameters:
   - N: the number of individuals
")set_random_genotype;



%ignore    update_traits;
%ignore    reset;
%ignore    point_sub_pop;
%ignore    set_outcrossing_rate;
%ignore    get_outcrossing_rate;
%ignore    set_crossover_rate;
%ignore    get_crossover_rate;
%ignore    set_recombination_model;
%ignore    get_recombination_model;
%ignore    update_fitness;
%ignore    point_sub_pop;

%pythoncode {
def random_genomes(self, n, location=-1):
    '''Get a sample of random genomes from the population

    Parameters:
       - n: number of random genomes to compute
       - location: location to sample genome from. default -1 leads to sampling from entire population

    Returns:
       - gts: (n x L) bool matrix with the n genotypes
    '''

    L = self._L()
    genotypes = _np.zeros((n, L), bool)
    for i in xrange(genotypes.shape[0]):
        genotypes[i] = self.get_genotype(self.get_random_clone(location))
    return genotypes
}

%feature("autodoc",
"Get a clonal genotype

Parameters: 
   - clone: clone 

Returns:
   - genotype: clonal genotype as numpy array of length L 
")get_genotype;

%feature("autodoc",
"Get random clone

Returns:
   - clone
")get_random_clone;
%exception get_random_clone {
        try {
                $action
        } catch (int err) {
            if (err == HP_BADARG)
                   PyErr_SetString(PyExc_ValueError,
                                "Location with this index does not exist!");
            else 
                   PyErr_SetString(PyExc_ValueError,
                                "Unknown exception is thrown. Aborting.");

            SWIG_fail;
        }
}
} /*extend multi_population*/ 
