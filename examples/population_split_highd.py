# vim: fdm=indent
'''
author:     Fabio Zanini, Richard Neher
date:       08/12/12
content:    Split a parent population into two daughter populations.
'''
# Import modules (setting the path should not be necessary when the module is
# installed in the PYTHONPATH)
import sys
sys.path.insert(0, '../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import FFPopSim as h


# specify parameters
N = 1000                            # Population size
L = 100                             # number of loci
mu = 0.000001                         # no new mutations
r = 0.01                            # crossover rate
outc = 0.01                         # outcrossing rate

### set up
pop = h.haploid_highd(L)            # produce an instance of haploid_lowd with L loci
pop.carrying_capacity = N           # set the steady-state population size
pop.outcrossing_rate = outc         # set the outcrossing rate
pop.crossover_rate = r              # assign the recombination rate
pop.mutation_rate = mu              # assign the mutation rate

# initialize the population with
#    - N/2 individuals with genotype 0, that is ----
#    - N/2 with the opposite genotype, that is ++++
pop.set_wildtype(N)

# Fitness
pop.set_trait_additive(0.001 * np.ones(L), 0)
pop.add_trait_coefficient(0.13, [0,2], 0)

# evolve for some time
pop.evolve(100)

# copy the population and evolve the two in parallel
pop2 = pop.copy()
g = []
af1 = []
af2 = []
for i in xrange(100):
    pop.evolve(10)
    pop2.evolve(10)
    af1.append(pop.get_allele_frequencies())
    af2.append(pop2.get_allele_frequencies())
    g.append(pop.generation)
af1 = np.array(af1)
af2 = np.array(af2)

# Plot the allele frequencies
colors = [cm.jet(int(255.0 * i / af1.shape[1])) for i in xrange(af1.shape[1])]
for i in xrange(L):
    if (af1[:, i] > 0).any() or (af2[:, i] > 0).any():
        plt.plot(g, af1[:, i], color=colors[i], ls='-')
        plt.plot(g, af2[:, i], color=colors[i], ls='--')

plt.xlabel('generations')
plt.ylabel('allele frequency')

plt.ion()
plt.show()

