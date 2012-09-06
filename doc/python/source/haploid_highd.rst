.. _haploid_highd:

.. currentmodule:: FFPopSim

haploid_highd Class Reference
=============================
.. autoclass:: haploid_highd
   :members: L,
             number_of_loci,
             N,
             population_size,
             generation,
             circular,
             carrying_capacity,
             recombination_model,
             outcrossing_rate,
             crossover_rate,
             mutation_rate,
             participation_ratio,
             number_of_clones,
             number_of_traits,
             status,
             max_fitness,
             set_allele_frequencies,
             set_genotypes,
             distance_Hamming,
             evolve,
             get_allele_frequencies,
             get_clone_sizes,
             get_fitnesses,
             get_trait_additive,
             get_genotype,
             get_genotypes,
             get_divergence_histogram,
             get_diversity_histogram,
             get_fitness_histogram,
             plot_divergence_histogram,
             plot_diversity_histogram,
             plot_fitness_histogram,
             random_genomes

   .. automethod:: haploid_highd.__init__
   .. automethod:: haploid_highd.add_trait_coefficient(value, loci, t=0)
   .. automethod:: haploid_highd.set_random_trait_epistasis(epistasis_std, t=0)
   .. automethod:: haploid_highd.add_fitness_coefficient(value, loci)
   .. automethod:: haploid_highd.set_random_epistasis(epistasis_std)
   .. automethod:: haploid_highd.add_genotype(genotype, n=1)
   .. automethod:: haploid_highd.set_fitness_additive(coefficients)
   .. automethod:: haploid_highd.set_trait_additive(coefficients, t)
   .. automethod:: haploid_highd.set_wildtype(N)
   .. automethod:: haploid_highd.flip_single_locus(locus)
   .. automethod:: haploid_highd.calc_stat()
   .. automethod:: haploid_highd.clear_trait(t=0)
   .. automethod:: haploid_highd.clear_fitness()
   .. automethod:: haploid_highd.clear_traits()
   .. automethod:: haploid_highd.bottleneck(size_of_bottleneck)
   .. automethod:: haploid_highd.get_allele_frequency(locus)
   .. automethod:: haploid_highd.get_chi(locus)
   .. automethod:: haploid_highd.get_chi2(locus1, locus2)
   .. automethod:: haploid_highd.get_LD(locus1, locus2)
   .. automethod:: haploid_highd.get_moment(locus1, locus2)
   .. automethod:: haploid_highd.get_pair_frequency(locus1, locus2)
   .. automethod:: haploid_highd.get_clone_size(n)
   .. automethod:: haploid_highd.get_divergence_statistics(n_sample=1000)
   .. automethod:: haploid_highd.get_diversity_statistics(n_sample=1000)
   .. automethod:: haploid_highd.get_fitness(n)
   .. automethod:: haploid_highd.get_fitness_statistics()
   .. automethod:: haploid_highd.get_trait(n, t=0)
   .. automethod:: haploid_highd.get_trait_statistics(t=0)
   .. automethod:: haploid_highd.get_trait_covariance(t1, t2)
   .. automethod:: haploid_highd.random_clone()
   .. automethod:: haploid_highd.random_clones(n)
   .. automethod:: haploid_highd.unique_clones()
