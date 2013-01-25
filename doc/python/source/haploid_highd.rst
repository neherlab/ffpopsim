.. _haploid_highd:

.. currentmodule:: FFPopSim

haploid_highd Class Reference
=============================
.. autoclass:: haploid_highd
   :members:

   .. automethod:: haploid_highd.__init__
   .. automethod:: haploid_highd.__str__
   .. automethod:: haploid_highd.__repr__
   .. automethod:: haploid_highd.copy
   .. automethod:: haploid_highd.status
   .. autoattribute:: haploid_highd.L
   .. autoattribute:: haploid_highd.number_of_loci
   .. autoattribute:: haploid_highd.N
   .. autoattribute:: haploid_highd.population_size
   .. autoattribute:: haploid_highd.generation
   .. autoattribute:: haploid_highd.number_of_clones
   .. autoattribute:: haploid_highd.number_of_traits
   .. autoattribute:: haploid_highd.max_fitness
   .. autoattribute:: haploid_highd.participation_ratio
   .. autoattribute:: haploid_highd.circular
   .. autoattribute:: haploid_highd.carrying_capacity
   .. autoattribute:: haploid_highd.recombination_model
   .. autoattribute:: haploid_highd.outcrossing_rate
   .. autoattribute:: haploid_highd.crossover_rate
   .. autoattribute:: haploid_highd.mutation_rate
   .. autoattribute:: haploid_highd.trait_weights
   .. automethod:: haploid_highd.get_clone(n)
   .. automethod:: haploid_highd.set_wildtype(N)
   .. automethod:: haploid_highd.set_allele_frequencies(frequencies, N)
   .. automethod:: haploid_highd.set_genotypes(genotypes, counts)
   .. automethod:: haploid_highd.evolve(gen=1)
   .. automethod:: haploid_highd.bottleneck(size_of_bottleneck)
   .. automethod:: haploid_highd.flip_single_locus(locus)
   .. automethod:: haploid_highd.calc_stat()
   .. automethod:: haploid_highd.get_divergence_statistics(n_sample)
   .. automethod:: haploid_highd.get_diversity_statistics(n_sample)
   .. automethod:: haploid_highd.get_trait_statistics(t)
   .. automethod:: haploid_highd.get_fitness_statistics()
   .. automethod:: haploid_highd.get_trait_covariance(t1, t2)
   .. automethod:: haploid_highd.get_allele_frequency(locus)
   .. automethod:: haploid_highd.get_allele_frequencies()
   .. automethod:: haploid_highd.get_pair_frequency(locus1, locus2)
   .. automethod:: haploid_highd.get_chi(locus)
   .. automethod:: haploid_highd.get_chi2(locus1, locus2)
   .. automethod:: haploid_highd.get_LD(locus1, locus2)
   .. automethod:: haploid_highd.get_moment(locus1, locus2)
   .. automethod:: haploid_highd.add_genotype(genotype, n=1)
   .. automethod:: haploid_highd.get_trait_additive(t)
   .. automethod:: haploid_highd.clear_trait(t)
   .. automethod:: haploid_highd.clear_traits()
   .. automethod:: haploid_highd.clear_fitness()
   .. automethod:: haploid_highd.set_trait_additive(coefficients, t)
   .. automethod:: haploid_highd.set_fitness_additive(coefficients)
   .. automethod:: haploid_highd.add_trait_coefficient(value, loci, t=0)
   .. automethod:: haploid_highd.add_fitness_coefficient(value, loci)
   .. automethod:: haploid_highd.get_trait_epistasis(t=0)
   .. automethod:: haploid_highd.set_random_trait_epistasis(epistasis_std, t=0)
   .. automethod:: haploid_highd.set_random_epistasis(epistasis_std)
   .. automethod:: haploid_highd.get_fitness(n)
   .. automethod:: haploid_highd.get_fitnesses
   .. automethod:: haploid_highd.get_trait(n, t=0)
   .. automethod:: haploid_highd.get_traits
   .. automethod:: haploid_highd.get_clone_size(n)
   .. automethod:: haploid_highd.get_clone_sizes
   .. automethod:: haploid_highd.get_genotype(n)
   .. automethod:: haploid_highd.get_genotypes
   .. automethod:: haploid_highd.unique_clones()
   .. automethod:: haploid_highd.distance_Hamming
   .. automethod:: haploid_highd.random_genomes
   .. automethod:: haploid_highd.random_clone()
   .. automethod:: haploid_highd.random_clones
   .. automethod:: haploid_highd.get_fitness_histogram
   .. automethod:: haploid_highd.plot_fitness_histogram
   .. automethod:: haploid_highd.get_divergence_histogram
   .. automethod:: haploid_highd.plot_divergence_histogram
   .. automethod:: haploid_highd.get_diversity_histogram
   .. automethod:: haploid_highd.plot_diversity_histogram
   .. automethod:: haploid_highd.genealogy
   .. automethod:: haploid_highd.track_locus_genealogy(loci)

