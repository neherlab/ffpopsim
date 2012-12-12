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
N = 500000                          # Population size
L = 4                               # number of loci
mu = 0.0                            # no new mutations
r = [0.01, 0.03, 0.05]              # recombination rates

### set up
pop = h.haploid_lowd(L)             # produce an instance of haploid_lowd with L loci
pop.carrying_capacity = N           # set the steady-state population size
pop.set_recombination_rates(r, h.CROSSOVERS)      # assign the recombination rate
pop.set_mutation_rates(mu)          # assign the mutation rate

# initialize the population with
#    - N/2 individuals with genotype 0, that is ----
#    - N/2 with the opposite genotype, that is ++++
pop.set_genotypes([0, 2**L - 1],
                  [N/2, N/2])

# evolve for some time
pop.evolve(50)

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
for i in xrange(af1.shape[1]):
    plt.plot(g, af1[:, i], color=colors[i], ls='-')
    plt.plot(g, af2[:, i], color=colors[i], ls='--')

plt.xlabel('generations')
plt.ylabel('allele frequency')

plt.ion()
plt.show()
