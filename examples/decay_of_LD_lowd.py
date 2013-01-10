# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       23/08/12
content:    Example of haploid_lowd on linkage relaxation via recombination
'''
# Import modules (setting the path should not be necessary when the module is
# installed in the PYTHONPATH)
import sys
sys.path.insert(0, '../pkg/python')

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
pop.set_recombination_rates(r)      # assign the recombination rate
pop.set_mutation_rates(mu)          # assign the mutation rate

# initialize the population with
#    - N/2 individuals with genotype 0, that is ----
#    - N/2 with the opposite genotype, that is ++++
pop.set_genotypes([0, 2**L - 1],
                  [N/2, N/2])

print "\nTrack LD and compare to deterministic expectations\n"

pop.status()

# get initial LD
LD_trajectories = [[pop.generation, pop.get_LD(0,1), pop.get_LD(0,2), pop.get_LD(0,3)]]

# evolve with accuracy of 5 generations and save LD along the way
for ii in xrange(50):

    pop.evolve(5)
    
    # get LD and time
    LD_trajectories.append([pop.generation, pop.get_LD(0,1), pop.get_LD(0,2), pop.get_LD(0,3)])

LD_trajectories=np.array(LD_trajectories)

# plot the LD trajectories and compare to exponential decay
cols = ['r', 'b', 'g', 'm', 'c']
for ii in xrange(LD_trajectories.shape[1]-1):
    plt.plot(LD_trajectories[:,0], LD_trajectories[:,ii+1], color=cols[ii], label=r'$D_{0'+str(ii+1)+'}$')
    plt.plot(LD_trajectories[:,0], 0.25 * np.exp(-LD_trajectories[:,0] * r * (ii+1)), ls='--', color=cols[ii])

plt.legend()
plt.xlabel('Time [generations]')
plt.ylabel('LD $D_{ij}$')

plt.ion()
plt.show()
