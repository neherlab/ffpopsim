.. _haploid_lowd:

.. currentmodule:: FFPopSim

haploid_lowd Class Reference
============================
.. autoclass:: haploid_lowd
   :members:

   .. automethod:: haploid_lowd.__init__
   .. automethod:: haploid_lowd.__str__
   .. automethod:: haploid_lowd.__repr__
   .. automethod:: haploid_lowd.copy
   .. automethod:: haploid_lowd.status
   .. autoattribute:: haploid_lowd.L
   .. autoattribute:: haploid_lowd.number_of_loci
   .. autoattribute:: haploid_lowd.N
   .. autoattribute:: haploid_lowd.population_size
   .. autoattribute:: haploid_lowd.generation
   .. autoattribute:: haploid_lowd.circular
   .. autoattribute:: haploid_lowd.carrying_capacity
   .. autoattribute:: haploid_lowd.recombination_model
   .. autoattribute:: haploid_lowd.outcrossing_rate
   .. automethod:: haploid_lowd.set_allele_frequencies(frequencies, N)
   .. automethod:: haploid_lowd.set_genotypes
   .. automethod:: haploid_lowd.set_wildtype(N)
   .. automethod:: haploid_lowd.set_recombination_rates
   .. automethod:: haploid_lowd.get_recombination_rates
   .. automethod:: haploid_lowd.set_mutation_rates
   .. automethod:: haploid_lowd.get_mutation_rates
   .. automethod:: haploid_lowd.evolve(gen=1)
   .. automethod:: haploid_lowd.evolve_deterministic(gen=1)
   .. automethod:: haploid_lowd.evolve_norec(gen=1)
   .. automethod:: haploid_lowd.get_genotype_frequency(genotype)
   .. automethod:: haploid_lowd.get_genotype_frequencies
   .. automethod:: haploid_lowd.get_allele_frequency(locus)
   .. automethod:: haploid_lowd.get_allele_frequencies
   .. automethod:: haploid_lowd.get_pair_frequency(locus1, locus2)
   .. automethod:: haploid_lowd.get_chi(locus)
   .. automethod:: haploid_lowd.get_chi2(locus1, locus2)
   .. automethod:: haploid_lowd.get_LD(locus1, locus2)
   .. automethod:: haploid_lowd.get_moment(locus1, locus2)
   .. automethod:: haploid_lowd.random_genomes
   .. automethod:: haploid_lowd.get_fitness(genotype)
   .. automethod:: haploid_lowd.get_fitnesses()
   .. automethod:: haploid_lowd.get_fitness_histogram
   .. automethod:: haploid_lowd.plot_fitness_histogram
   .. automethod:: haploid_lowd.get_divergence_statistics
   .. automethod:: haploid_lowd.get_divergence_histogram
   .. automethod:: haploid_lowd.plot_divergence_histogram
   .. automethod:: haploid_lowd.get_diversity_statistics
   .. automethod:: haploid_lowd.get_diversity_histogram
   .. automethod:: haploid_lowd.plot_diversity_histogram
   .. automethod:: haploid_lowd.set_fitness_function
   .. automethod:: haploid_lowd.set_fitness_coefficients
   .. automethod:: haploid_lowd.set_fitness_additive(coefficients)
   .. automethod:: haploid_lowd.genotype_entropy()
   .. automethod:: haploid_lowd.allele_entropy()


