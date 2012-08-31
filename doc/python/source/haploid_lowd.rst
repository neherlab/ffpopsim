.. _haploid_lowd:

.. currentmodule:: FFPopSim

haploid_lowd Class Reference
============================
.. autoclass:: haploid_lowd
   :members: L,
             number_of_loci,
             N,
             population_size,
             generation,
             circular,
             carrying_capacity,
             recombination_model,
             outcrossing_rate,
             random_genomes,
             set_allele_frequencies,
             set_fitness_function,
             set_genotypes,
             set_recombination_rates,
             set_mutation_rates,
             get_allele_frequencies,
             get_fitnesses,
             get_genotype_frequencies,
             get_mutation_rates,
             get_diversity_statistics,
             get_diversity_histogram,
             get_divergence_statistics,
             get_divergence_histogram,
             get_fitness_statistics,
             get_fitness_histogram,
             plot_diversity_histogram,
             plot_divergence_histogram,
             plot_fitness_histogram

   .. automethod:: haploid_lowd.__init__
   .. automethod:: haploid_lowd.allele_entropy()
   .. automethod:: haploid_lowd.evolve(gen=1)
   .. automethod:: haploid_lowd.evolve_deterministic(gen=1)
   .. automethod:: haploid_lowd.evolve_norec(gen=1)
   .. automethod:: haploid_lowd.genotype_entropy()
   .. automethod:: haploid_lowd.get_allele_frequency(locus)
   .. automethod:: haploid_lowd.get_pair_frequency(locus1, locus2)
   .. automethod:: haploid_lowd.get_chi(locus)
   .. automethod:: haploid_lowd.get_chi2(locus1, locus2)
   .. automethod:: haploid_lowd.get_fitness(genotype)
   .. automethod:: haploid_lowd.get_genotype_frequency(genotype)
   .. automethod:: haploid_lowd.get_LD(locus1, locus2)
   .. automethod:: haploid_lowd.get_moment(locus1, locus2)
   .. automethod:: haploid_lowd.set_fitness_additive(coefficients)
   .. automethod:: haploid_lowd.set_wildtype(N)


