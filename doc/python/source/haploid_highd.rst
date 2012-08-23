.. _haploid_highd:

haploid_highd Class Reference
=============================
.. autoclass:: FFPopSim.haploid_highd
   :members: L,
             N,
             number_of_loci,
             population_size,
             circular,
             carrying_capacity,
             mutation_rate,
             outcrossing_rate,
             crossover_rate,
             recombination_model,
             generation,
             distance_Hamming,
             evolve,
             get_additive_trait,
             get_allele_frequencies,
             get_clone_sizes,
             get_divergence_histogram,
             get_diversity_histogram,
             get_fitness_histogram,
             plot_divergence_histogram,
             plot_diversity_histogram,
             plot_fitness_histogram,
             get_fitnesses, get_genotype,
             get_genotypes,
             random_genomes,
             participation_ratio,
             number_of_clones,
             number_of_traits,
             set_allele_frequencies,
             set_genotypes,
             max_fitness

   .. automethod:: FFPopSim.haploid_highd.add_trait_coefficient(value, loci, t)
   .. automethod:: FFPopSim.haploid_highd.add_fitness_coefficient(value, loci)
   .. automethod:: FFPopSim.haploid_highd.add_genotypes(gt, n=1)
   .. automethod:: FFPopSim.haploid_highd.bottleneck(size_of_bottleneck)
   .. automethod:: FFPopSim.haploid_highd.calc_stat()
   .. automethod:: FFPopSim.haploid_highd.clear_trait(t)
   .. automethod:: FFPopSim.haploid_highd.clear_fitness()
   .. automethod:: FFPopSim.haploid_highd.clear_traits()
   .. automethod:: FFPopSim.haploid_highd.get_allele_frequency(locus)
   .. automethod:: FFPopSim.haploid_highd.get_chi(locus)
   .. automethod:: FFPopSim.haploid_highd.get_pair_frequency(locus1, locus2)
   .. automethod:: FFPopSim.haploid_highd.get_clone_size(n)
   .. automethod:: FFPopSim.haploid_highd.get_divergence_statistics(n_sample=1000)
   .. automethod:: FFPopSim.haploid_highd.get_diversity_statistics(n_sample=1000)
   .. automethod:: FFPopSim.haploid_highd.get_fitness(n)
   .. automethod:: FFPopSim.haploid_highd.get_fitness_statistics()
   .. automethod:: FFPopSim.haploid_highd.get_trait(n, t=0)
   .. automethod:: FFPopSim.haploid_highd.get_trait_statistics(t=0)
   .. automethod:: FFPopSim.haploid_highd.get_trait_covariance(t1, t2)
   .. automethod:: FFPopSim.haploid_highd.random_clone()
   .. automethod:: FFPopSim.haploid_highd.random_clones(n)
   .. automethod:: FFPopSim.haploid_highd.set_additive_fitness(coefficients)
   .. automethod:: FFPopSim.haploid_highd.set_additive_trait(coefficients, t)
   .. automethod:: FFPopSim.haploid_highd.set_random_epistasis(epistasis_std)
   .. automethod:: FFPopSim.haploid_highd.set_random_trait_epistasis(epistasis_std)
   .. automethod:: FFPopSim.haploid_highd.set_wildtype(N)
   .. automethod:: FFPopSim.haploid_highd.unique_clones()
