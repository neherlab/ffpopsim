# vim: fdm=indent
'''
author:     Richard Neher
date:       23/08/12
content:    Example of haploid_lowd on linkage relaxation via recombination
'''
# Import module
import sys
sys.path.append('../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h


# specify parameters
N = 500000                          # Population size
L = 4                               # number of loci
mu = 0.0                            # no new mutations
r = 0.01                            # recombination rate

### set up
pop = h.haploid_lowd(L)             # produce an instance of haploid_lowd with L loci
pop.carrying_capacity = N           # set the steady-state population size
pop.set_recombination_rates(r)      # assign the recombination rates
pop.set_mutation_rates(mu)          # assign the mutation rate

# initialize the population with N/2 individuals with genotypes 0, that is ----
# and N/2 with the opposite genotype, that is ++++
pop.set_genotypes([0, 2**L-1],[N/2, N/2])


max_gen = 50
LD_trajectories = [[pop.generation,pop.get_LD(0,1), pop.get_LD(0,2), pop.get_LD(0,3)]]
for ii in xrange(max_gen):
    pop.evolve(5)               #N/10 generations between successive samples
    LD_trajectories.append([pop.generation, pop.get_LD(0,1), pop.get_LD(0,2), pop.get_LD(0,3)])

LD_trajectories=np.array(LD_trajectories)
cols = ['r', 'b', 'g', 'm', 'c']
for ii in xrange(LD_trajectories.shape[1]-1):
    plt.plot(LD_trajectories[:,0], LD_trajectories[:,ii+1], color=cols[ii], label=r'$D_{0'+str(ii+1)+'}$')
    plt.plot(LD_trajectories[:,0], 0.25 * np.exp(-LD_trajectories[:,0]* r * (ii+1)), ls='--', color=cols[ii])

plt.legend()
plt.xlabel('Time [generations]')
plt.ylabel('LD $D_{ij}$')

plt.ion()
plt.show()
