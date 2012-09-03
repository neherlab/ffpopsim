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
N = 1e12                        # Population size
Lmax_general = 13               # Maximal number of loci
Lmax_single_xo = 17             # Maximal number of loci
r = 0.01                        # Recombination rate
mu = 0.001                      # Mutation rate
G = 1000                         # Generations

# Repeat the same simulation for various numbers of loci, and see how the
# algorithm scales. It should be O(3^L) with recombination, and O(L 2^L) with
# single crossovers only or skipping recombination altogether.
print "general"
exec_time = []
for L in range(2,Lmax_general+1):
    t1=time.time()
    ### set up
    pop = h.haploid_lowd(L)     # produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = N   # set the population size

    # set and additive fitness function with random coefficients. 
    # Note that FFPopSim models fitness landscapes in the -1/+1 basis
    pop.set_fitness_additive(0.01*np.random.randn(L))

    pop.set_recombination_rates(r)  # assign the recombination rates
    pop.set_mutation_rates(mu)  # assign the mutation rate
    
    #initialize the population with N individuals in linkage equilibrium
    pop.set_allele_frequencies(0.2*np.ones(L), N)

    pop.evolve(G)                   # run for G generations to measure execution time
    
    t2=time.time()

    exec_time.append([L, t2-t1])    # store the execution time
    
exec_time=np.array(exec_time)

# repeat with single crossover recombination
print "single crossover"
exec_time_single_xo = []
for L in range(2,Lmax_single_xo+1):
    t1=time.time()
    ### set up
    pop = h.haploid_lowd(L)     # produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = N   # set the population size
    pop.set_fitness_additive(0.01*np.random.randn(L))
    # assign the recombination rates, assume SINGLE CROSSOVER (otherwise everything is the same)
    pop.set_recombination_rates(r, h.SINGLE_CROSSOVER)  
    pop.set_mutation_rates(mu)  # assign the mutation rate
    pop.set_allele_frequencies(0.2*np.ones(L), N)

    pop.evolve(G)               # run for G generations to measure execution time
    
    t2=time.time()

    exec_time_single_xo.append([L, t2-t1])    # store the execution time
    
exec_time_single_xo=np.array(exec_time_single_xo)


# repeat without recombination
print "no recombination"
exec_time_norec = []
for L in range(2,Lmax_single_xo+1):
    t1=time.time()
    ### set up
    pop = h.haploid_lowd(L)
    pop.carrying_capacity = N
    pop.set_fitness_additive(0.01*np.random.randn(L)) 
    pop.set_mutation_rates(mu)
    pop.set_allele_frequencies(0.2*np.ones(L), N)
    pop.evolve_norec(G)         #evolve for G generation skipping the recombination step
    t2=time.time()
    exec_time_norec.append([L, t2-t1])

exec_time_norec=np.array(exec_time_norec)

# Plot the execution times and lines indicating the expected computational complexity
plt.figure()
#general recombination
plt.plot(exec_time[:,0], exec_time[:,1],label='with recombination', linestyle='None', marker = 'o')
plt.plot(exec_time[:,0], exec_time[-1,1]/3.0**(Lmax_general-exec_time[:,0]),label=r'$\propto 3^L$')

##single crossovers
plt.plot(exec_time_single_xo[:,0],
         exec_time_single_xo[:,1],label='single crossover', linestyle='None', marker = 'o')
plt.plot(exec_time_single_xo[:,0],
         exec_time_single_xo[-1,1]/2.0**(Lmax_single_xo-exec_time_single_xo[:,0])*(exec_time_single_xo[:,0]/Lmax_single_xo),
         label=r'$\propto L2^L$')

# no recombination
plt.plot(exec_time_norec[:,0], exec_time_norec[:,1],label='without recombination', linestyle='None', marker = 'x')

# complexity of the naive algorithm
plt.plot(exec_time[:,0], exec_time[-1,1]/3.0**(Lmax_general)*8**(exec_time[:,0]),label=r'$\propto 8^L$')

ax=plt.gca()
ax.set_yscale('log')
plt.xlabel('number of loci')
plt.ylabel('seconds for '+str(G)+' generations')
plt.legend(loc=4)
plt.xlim([1,Lmax_single_xo+2])
plt.ylim([0.5*np.min(exec_time_norec[:,1]),3*np.max(exec_time[:,1])])


plt.ion()
plt.show()
