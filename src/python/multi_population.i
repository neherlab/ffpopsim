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

%rename(_get_migration_rate) get_migration_rate;
%rename(_set_migration_rate) set_migration_rate;
%pythoncode{
@property
def migration_rate(self):
   '''Populational migration rate'''
   return self._get_migration_rate()

@migration_rate.setter
def migration_rate(self, m):
        self._set_migration_rate(m)
}


%rename(_get_number_of_locations) get_locations;
%pythoncode{
@property
def locations_number(self):
    '''Number of locations (read-only)'''
    return self._get_number_of_locations()
}


%rename(_max_fitness) max_fitness;
%pythoncode{
@property
def max_fitness(self):
    '''Maximal fitness of the population (read-only)'''
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

%rename (_set_carrying_capacity) set_carrying_capacity;
%rename (_get_carrying_capacity) get_carrying_capacity;
%pythoncode{
@property
def carrying_capacity(self):
    '''Carrying capacity of the population'''
    return self._get_carrying_capacity()
@carrying_capacity.setter
def carrying_capacity(self, N):
    self._set_carrying_capacity(N)
}

%rename(_N) N;
%pythoncode{
@property
def N(self):
    '''population size (read-only)'''
    return self._N()
}

%rename(_set_mutation_rate) set_mutation_rate;
%rename(_get_mutation_rate) get_mutation_rate;
%pythoncode
{
@property
def mutation_rate(self):
    '''Mutation rate per locus'''
    return self._get_mutation_rate()
@mutation_rate.setter
def mutation_rate(self, mu):
    self._set_mutation_rate(mu)
}

%ignore set_trait_coefficients;
%ignore set_trait_weights;


/* genealogy */
%feature("autodoc",
"Track the genealogy of some loci.

   Parameters:
   - loci: sites whose genealogy is being stored

   Returns:
   - zero if successful, error code otherwise
") track_locus_genealogy;


%ignore    genealogy;


/*evolution*/
%feature("autodoc",
"Evolve for some generations.

Parameters:
   - gen: number of generations, defaults to one
") evolve;



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
%ignore   point_sub_pop;

} /*extend multi_population*/ 
