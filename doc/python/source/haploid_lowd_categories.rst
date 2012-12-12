.. currentmodule:: FFPopSim

haploid_lowd Overview
=====================
The following lists all methods and attributes of ``haploid_lowd``, grouped according to purpose and function.
Member functions are shown with the prefix ''haploid_lowd.*''. If you have initialized for example as
``pop=haploid_lowd(L)``, the prefix in your program needs to be ``pop.*``.

.. note:: Clicking on the name of the function will take you to a more detailed explanation listing all arguments.

Initialization
^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.__init__
   haploid_lowd.copy

Attributes
^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.L
   haploid_lowd.number_of_loci
   haploid_lowd.N
   haploid_lowd.population_size
   haploid_lowd.generation

   haploid_lowd.carrying_capacity
   haploid_lowd.circular
   haploid_lowd.outcrossing_rate
   haploid_lowd.recombination_model

Status
^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.status

Initialize the Population
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.set_wildtype
   haploid_lowd.set_allele_frequencies
   haploid_lowd.set_genotypes

Set the fitness landscape
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.set_fitness_additive
   haploid_lowd.set_fitness_function
   haploid_lowd.set_fitness_coefficients


Mutation and Recombination
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_lowd.circular
   haploid_lowd.outcrossing_rate
   haploid_lowd.recombination_model

   haploid_lowd.set_recombination_rates
   haploid_lowd.get_recombination_rates
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

Analyze the fitness distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:


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


