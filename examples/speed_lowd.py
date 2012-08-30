# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       24/08/12
content:    Example of haploid_lowd on the algorithm complexity
'''
# Import module
import sys
sys.path.append('../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h
import time

# specify parameters
N = 1000                        # Population size
Lmax = 15                       # Maimal number of loci
r = 0.01                        # Recombination rate
mu = 0.001                      # Mutation rate
G = 1000                        # Generations

# Repeat the same simulation for various numbers of loci, and see how the
# algorithm scales. It should be O(3^L) with recombination, and O(L 2^L) without
# recombination.
exec_time = []
for L in range(2,Lmax-2):
    t1=time.time()
    ### set up
    pop = h.haploid_lowd(L)     # produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = N   # set the population size

    # set and additive fitness function. Note that FFPopSim models fitness landscape
    pop.set_fitness_additive(0.01*np.random.randn(L))

    pop.set_recombination_rates(r)  # assign the recombination rates
    pop.set_mutation_rates(mu)  # assign the mutation rate
    
    #initialize the population with N individuals with genotypes 0, that is ----
    pop.set_allele_frequencies(0.2*np.ones(L), N)

    pop.evolve(G)               # run for G generations to equilibrate
    
    t2=time.time()

    exec_time.append([L, t2-t1])    # store the execution time
    
exec_time=np.array(exec_time)

# repeat with single crossover recombination
exec_time_sc = []
for L in range(2,Lmax):
    t1=time.time()
    ### set up
    pop = h.haploid_lowd(L)     # produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = N   # set the population size

    # set and additive fitness function. Note that FFPopSim models fitness landscape
    pop.set_fitness_additive(0.01*np.random.randn(L))

    pop.set_recombination_rates(r, h.SINGLE_CROSSOVER)  # assign the recombination rates
    pop.set_mutation_rates(mu)  # assign the mutation rate
    
    #initialize the population with N individuals with genotypes 0, that is ----
    pop.set_allele_frequencies(0.2*np.ones(L), N)

    pop.evolve(G)               # run for G generations to equilibrate
    
    t2=time.time()

    exec_time_sc.append([L, t2-t1])    # store the execution time
    
exec_time_sc=np.array(exec_time_sc)


# repeat without recombination
exec_time_norec = []
for L in range(2,Lmax):
    t1=time.time()
    ### set up
    pop = h.haploid_lowd(L)
    pop.carrying_capacity = N
    pop.set_fitness_additive(0.01*np.random.randn(L)) 
    pop.set_mutation_rates(mu)
    pop.set_allele_frequencies(0.2*np.ones(L), N)
    pop.evolve_norec(G)
    t2=time.time()
    exec_time_norec.append([L, t2-t1])

exec_time_norec=np.array(exec_time_norec)

# Plot the execution times
plt.figure()
plt.plot(exec_time[:,0], exec_time[:,1],label='with recombination', linestyle='None', marker = 'o')
plt.plot(exec_time[:,0], exec_time[-1,1]/3.0**(Lmax-2-exec_time[:,0]-1),label=r'$\propto 3^L$')

plt.plot(exec_time_sc[:,0],
         exec_time_sc[:,1],label='single crossover', linestyle='None', marker = 'o')
plt.plot(exec_time_sc[:,0],
         exec_time_sc[-1,1]/2.0**(Lmax-exec_time_sc[:,0]-1),
         label=r'$\propto 2^L$')

plt.plot(exec_time_norec[:,0], exec_time_norec[:,1],label='without recombination', linestyle='None', marker = 'x')
plt.plot(exec_time_norec[:,0], exec_time_norec[-1,1]/2.0**(Lmax-exec_time_norec[:,0]-1),label=r'$\propto 2^L$')

plt.plot(exec_time[:,0], exec_time[-1,1]/3.0**(Lmax)*8**(exec_time[:,0]),label=r'$\propto 8^L$')

ax=plt.gca()
ax.set_yscale('log')
plt.xlabel('number of loci')
plt.ylabel('seconds for '+str(G)+' generations')
plt.legend(loc=2)
plt.xlim([1,Lmax])
plt.ylim([0.2*np.min(exec_time_norec[:,1]),10*np.max(exec_time[:,1])])

plt.ion()
plt.show()
