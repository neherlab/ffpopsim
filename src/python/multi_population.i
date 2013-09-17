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
  
/*Interface file for multi_population.h */

/*FIX FIXME actual class documentation here */
%define DOCSTRING_MULTI_POPULATION
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
   pop.mutation_rate = 1e-4       # mutation rate per site per generation
   pop.outcrossing_rate = 1e-1    # probability of sexual reproduction per gen
   pop.crossover_rate = 1e-2      # probability of crossover per site per gen
   pop.evolve(100)                # evolve for 100 generations
   c.plot_divergence_histogram()
   plt.show()
   ######################################

Populations can have a number of phenotypic traits that contribute to the fitness
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
%feature("autodoc",DOCSTRING_MULTI_POPULATION) multi_population;

%extend multi_population {
/*constructor*/
%feature("autodoc",
"Construct a high-dimensional structured population with certain parameters.
Parameters:
   - new_locations: number of locations 
   - L: number of loci
   - rng_seed: seed for the random generator. If zero (default) pick a random number
   - number_of_traits: number of phenotypic traits, defaults to one
   - all_polymorphic: option to use an infinite-sites model tracking ancestral alleles
                       (only available with a single phenotypic trait and zero mutation rate)
") multi_population;


%ignore    reset;
%ignore    point_sub_pop;
%ignore    get_locations;
%ignore    N;
%ignore    get_generation;
%ignore    max_fitness;
%ignore    number_of_migration_events;
%ignore    set_migration_rate;
%ignore    get_migration_rate;
%ignore    set_mutation_rate;
%ignore    get_mutation_rate;
%ignore    set_carrying_capacity;
%ignore    get_carrying_capacity;
%ignore    set_outcrossing_rate;
%ignore    get_outcrossing_rate;
%ignore    set_crossover_rate;
%ignore    get_crossover_rate;
%ignore    set_recombination_model;
%ignore    get_recombination_model;
%ignore    set_trait_coefficient;
%ignore    set_trait_weights;
%ignore    track_locus_genealogy;
%ignore    print_newick;
%ignore    evolve;
%ignore    migrate;
%ignore    migrate;
%ignore    add_random_genotype;
}/*extend multi_population*/

