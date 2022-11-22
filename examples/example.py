import FFPopSim as h

c = h.haploid_lowd(4)
c.set_allele_frequencies([0,0.3,0.6,0.9], N=1000)
c.evolve(100)
c.plot_diversity_histogram()
