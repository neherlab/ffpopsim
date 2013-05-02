.. currentmodule:: FFPopSim

haploid_highd Overview
=====================
The following lists all methods and attributes of ``haploid_highd``, grouped according to purpose and function.
Member functions are shown with the prefix ``haploid_highd.*``.
If you have initialized for example as ``pop=haploid_highd(L)``, the prefix in your program needs to be ``pop.*``.

.. note:: Clicking on the name of the function will take you to a more detailed explanation listing all arguments.

Initialization, copy, and storage
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.__init__
   haploid_highd.copy
   haploid_highd.dump

   load_haploid_highd


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
   haploid_highd.carrying_capacity
   haploid_highd.circular
   haploid_highd.outcrossing_rate
   haploid_highd.crossover_rate
   haploid_highd.recombination_model
   haploid_highd.mutation_rate
   haploid_highd.trait_weights

Status
^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.status

Initialize the population
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.set_wildtype
   haploid_highd.set_allele_frequencies
   haploid_highd.set_genotypes
   haploid_highd.add_genotype

Set phenotypes and fitness function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.trait_weights
   haploid_highd.set_trait_additive
   haploid_highd.add_trait_coefficient
   haploid_highd.set_random_trait_epistasis
   haploid_highd.clear_trait
   haploid_highd.clear_traits

The following methods are shortcuts if you have only one trait and will not work
for populations with more than one trait:

.. autosummary::
   :nosignatures:

   haploid_highd.set_fitness_additive
   haploid_highd.add_fitness_coefficient
   haploid_highd.set_random_epistasis
   haploid_highd.clear_fitness

Evolution
^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.evolve
   haploid_highd.bottleneck
   haploid_highd.unique_clones
   haploid_highd.flip_single_locus
   haploid_highd.calc_stat


Random Sampling
^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.random_clone
   haploid_highd.random_clones
   haploid_highd.random_genomes

Clone structure and genotypes 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.get_clone_size
   haploid_highd.get_clone_sizes
   haploid_highd.get_genotype
   haploid_highd.get_genotypes


Get allele frequencies and linkage disequilibria
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.get_allele_frequency
   haploid_highd.get_allele_frequencies
   haploid_highd.get_pair_frequency
   haploid_highd.get_chi
   haploid_highd.get_LD
   haploid_highd.get_chi2
   haploid_highd.get_moment

Analyze fitness and trait distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.calc_stat
   haploid_highd.get_trait_additive
   haploid_highd.get_trait_epistasis
   haploid_highd.get_trait_statistics
   haploid_highd.get_trait_covariance
   haploid_highd.get_fitness
   haploid_highd.get_fitnesses
   haploid_highd.get_fitness_statistics
   haploid_highd.get_fitness_histogram
   haploid_highd.plot_fitness_histogram


Divergence and Diversity
^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:
   
   haploid_highd.get_divergence_statistics
   haploid_highd.get_diversity_statistics
   haploid_highd.get_divergence_histogram
   haploid_highd.get_diversity_histogram
   haploid_highd.plot_divergence_histogram
   haploid_highd.plot_diversity_histogram


Genealogies and Trees
^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:

   haploid_highd.track_locus_genealogy
   haploid_highd.genealogy
