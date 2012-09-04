# vim: fdm=indent
'''
author:     Richard Neher
date:       23/08/12
content:    Example of haploid_lowd on linkage relaxation via recombination
'''
# Import module
import sys
sys.path.insert(0,'../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h


# specify parameters
N = 10000                         # Population size
L = 100                             # number of loci
mu = 0.0                            # no new mutations
r = 0.1/L                           # crossover rate

### set up
pop = h.haploid_highd(L)             # produce an instance of haploid_lowd with L loci
pop.carrying_capacity = N            # set the steady-state population size
pop.outcrossing_rate=1               # assign the outcrossing rate (equals 1 be default)
pop.crossover_rate = r               # assign the crossover rate per locus
pop.mutation_rate=mu                 # assign the mutation rate per locus

# initialize the population with N/2 individuals with genotypes 00000000...0000
# and N/2 with the opposite genotype 111111111111....111
pop.set_genotypes([np.zeros(L,dtype='int'), np.ones(L,dtype='int')],[N/2, N/2])

#locus pairs for which LD is to be tracked
locus_pairs = [ [0,10], [0,20], [40,50], [10,60]]

max_gen = 50
#get initial LD
LD_trajectories = [[pop.get_LD(l1,l2) for l1,l2 in locus_pairs]]
tp = [pop.generation]
for ii in xrange(max_gen):
    pop.evolve(5)               #5 generations between successive samples
    #get LD and time
    LD_trajectories.append([pop.get_LD(l1,l2) for l1,l2 in locus_pairs])
    tp.append(pop.generation)
    
LD_trajectories=np.array(LD_trajectories); tp=np.array(tp)

#plot the LD trajectories and compare to exponential decay
cols = ['r', 'b', 'g', 'm', 'c']
plt.figure()
for ii,(l1,l2) in enumerate(locus_pairs):
    plt.plot(tp,LD_trajectories[:,ii], color=cols[ii], label=r'$D_{'+str(l1)+','+str(l2)+'}$')
    plt.plot(tp,0.25*np.exp(-tp* r * abs(l1-l2)), ls='--', color=cols[ii])

plt.legend()
plt.title('Decay of LD and comparison to theory (dashed lines)')
plt.xlabel('Time [generations]')
plt.ylabel('LD $D_{ij}$')

plt.ion()
plt.show()

