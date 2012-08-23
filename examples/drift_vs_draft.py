import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h

L=256

### set up
pop = h.haploid_highd(L)             #produce an instance of haploid_lowd with L loci
pop.carrying_capacity = 10000
pop.outcrossing_rate = 0.1
pop.crossover_rate = 1.0/pop.L
pop.mutation_rate = 0.4/pop.carrying_capacity

selection_coefficients = 0.0*np.ones(pop.L)
selection_coefficients[::10] = -0.1

initial_allele_frequencies = 0.5*np.ones(pop.L)
initial_allele_frequencies[::10] = 0

pop.set_allele_frequencies(initial_allele_frequencies, pop.carrying_capacity)
pop.set_additive_trait(selection_coefficients)


maxgen = 1000
allele_frequencies = [pop.get_allele_frequencies()]
tp = [pop.generation]
while pop.generation<maxgen:
    pop.evolve(10)
    print pop.N
    if (pop.generation%50==0):
        selection_coefficients[10*np.random.randint(0,25)] = 0.01
        pop.set_additive_trait(selection_coefficients)
    allele_frequencies.append(pop.get_allele_frequencies())
    tp.append(pop.generation)


allele_frequencies=np.array(allele_frequencies)

plt.plot(tp, allele_frequencies)

