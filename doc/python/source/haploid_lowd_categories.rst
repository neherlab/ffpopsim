haploid_lowd Overview
=====================
``haploid_lowd`` has several kinds of methods.

.. currentmodule:: FFPopSim

Initialization
^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.__init__

Attributes
^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.L
   haploid_lowd.number_of_loci
   haploid_lowd.N
   haploid_lowd.population_size
   haploid_lowd.generation

   haploid_lowd.circular
   haploid_lowd.carrying_capacity
   haploid_lowd.free_recombination
   haploid_lowd.outcrossing_rate

Set Population
^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.set_wildtype
   haploid_lowd.set_allele_frequencies
   haploid_lowd.set_genotypes

Mutation and Recombination
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.set_recombination_rates
   haploid_lowd.set_mutation_rates
   haploid_lowd.get_mutation_rates

Evolution
^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.evolve
   haploid_lowd.evolve_deterministic
   haploid_lowd.evolve_norec

Random Sampling
^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.random_genomes

Get allele/genotype frequencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.get_allele_frequency
   haploid_lowd.get_allele_frequencies
   haploid_lowd.get_genotype_frequency
   haploid_lowd.get_genotype_frequencies
   haploid_lowd.get_pair_frequency
   haploid_lowd.get_chi
   haploid_lowd.get_chi2
   haploid_lowd.get_LD
   haploid_lowd.get_moment

   haploid_lowd.allele_entropy
   haploid_lowd.genotype_entropy

Fitness Landscape
^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.set_fitness_additive
   haploid_lowd.set_fitness_function
   haploid_lowd.get_fitness
   haploid_lowd.get_fitnesses

   haploid_lowd.get_fitness_statistics
   haploid_lowd.get_fitness_histogram
   haploid_lowd.plot_fitness_histogram

Divergence and Diversity
^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.get_divergence_statistics
   haploid_lowd.get_diversity_statistics

   haploid_lowd.get_divergence_histogram
   haploid_lowd.get_diversity_histogram

   haploid_lowd.plot_divergence_histogram
   haploid_lowd.plot_diversity_histogram


