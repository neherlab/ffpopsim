haploid_highd Overview
=====================
.. currentmodule:: FFPopSim

``haploid_highd`` has several kinds of methods.

Initialization
^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.__init__

Attributes
^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.L
   haploid_highd.number_of_loci
   haploid_highd.N
   haploid_highd.population_size
   haploid_highd.number_of_traits
   haploid_highd.generation
   haploid_highd.participation_ratio
   haploid_highd.number_of_clones

   haploid_highd.circular
   haploid_highd.carrying_capacity
   haploid_highd.recombination_model
   haploid_highd.outcrossing_rate
   haploid_highd.crossover_rate
   haploid_highd.mutation_rate

Set Population
^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.set_wildtype
   haploid_highd.set_allele_frequencies
   haploid_highd.set_genotypes
   haploid_highd.add_genotypes

Evolution
^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.evolve
   haploid_highd.bottleneck

Random Sampling
^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.random_clone
   haploid_highd.random_clones
   haploid_highd.random_genomes

Clone Structure
^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.get_clone_size
   haploid_highd.get_clone_sizes
   haploid_highd.get_genotype
   haploid_highd.get_genotypes
   haploid_highd.unique_clones

Get allele/genotype frequencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.get_allele_frequency
   haploid_highd.get_allele_frequencies
   haploid_highd.get_pair_frequency
   haploid_highd.get_chi
   haploid_highd.get_chi2
   haploid_highd.get_LD
   haploid_highd.get_moment

Fitness and Trait Landscapes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.get_fitness
   haploid_highd.get_fitnesses

   haploid_highd.set_trait_additive
   haploid_highd.get_trait_additive
   haploid_highd.add_trait_coefficient
   haploid_highd.set_random_trait_epistasis
   haploid_highd.get_trait_covariance
   haploid_highd.clear_trait
   haploid_highd.clear_traits

   haploid_highd.get_fitness_statistics
   haploid_highd.get_fitness_histogram
   haploid_highd.plot_fitness_histogram

The following methods are shortcuts if you have only one trait
(true if you do not subclass ``haploid_highd``) and will not work
for populations with more than one phenotypic trait:

.. autosummary::
   :nosignatures:

   haploid_highd.set_fitness_additive
   haploid_highd.add_fitness_coefficient
   haploid_highd.set_random_epistasis
   haploid_highd.clear_fitness

Divergence and Diversity
^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:
   
   haploid_highd.calc_stat

   haploid_highd.get_divergence_statistics
   haploid_highd.get_diversity_statistics
   haploid_highd.get_trait_statistics

   haploid_highd.get_divergence_histogram
   haploid_highd.get_diversity_histogram

   haploid_highd.plot_divergence_histogram
   haploid_highd.plot_diversity_histogram
