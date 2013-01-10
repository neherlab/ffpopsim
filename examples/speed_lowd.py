# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       24/08/12
content:    Example of haploid_lowd on the algorithm runtime complexity.
            The same series of simulations are performed:
                - with multiple crossovers, O(3^L),
                - with at most one crossover, O(L 2^L),
                - without recombination, O(L 2^L),
            in this order.
'''
# Import modules (setting the path should not be necessary when the module is
# installed in the PYTHONPATH)
import sys
sys.path.insert(0, '../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h
import time


# specify parameters
N = 1e12                        # population size
Lmax_general = 13               # maximal number of loci
Lmax_single_xo = 17             # maximal number of loci
r = 0.01                        # recombination rate
mu = 0.001                      # mutation rate
G = 100                         # generations

# 1. multiple crossovers
print "\ngeneral recombination"
exec_time = []
for L in range(2,Lmax_general+1):
    print L,"loci, maximum",Lmax_general,
    # set up population
    pop = h.haploid_lowd(L)     # produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = N   # set the population size

    # set and additive fitness function with random coefficients. 
    # Note that FFPopSim models fitness landscapes in the -1/+1 basis
    pop.set_fitness_additive(0.01 * np.random.randn(L))

    pop.set_recombination_rates(r)  # recombination rates
    pop.set_mutation_rates(mu)      # mutation rate
    
    # initialize the population with N individuals in linkage equilibrium
    pop.set_allele_frequencies(0.2 * np.ones(L), N)

    # run for G generations to measure execution time
    t1=time.time()
    pop.evolve(G)
    t2=time.time()

    print "time required for",G,"generations:",round(t2-t1,3),'s'
    exec_time.append([L, t2-t1])            # store the execution time
exec_time=np.array(exec_time)


# 2. single crossover recombination
print "\nsingle crossover"
exec_time_single_xo = []
for L in range(2,Lmax_single_xo+1):
    print L,"loci, maximum",Lmax_single_xo,

    # set up population
    pop = h.haploid_lowd(L)     # produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = N   # set the population size

    # set and additive fitness function with random coefficients. 
    pop.set_fitness_additive(0.01 * np.random.randn(L))
    
    # assign the recombination rates, assume SINGLE CROSSOVER (otherwise everything is the same)
    pop.set_recombination_rates(r, h.SINGLE_CROSSOVER)  
    pop.set_mutation_rates(mu)      # mutation rate

    # initialize the population with N individuals in linkage equilibrium
    pop.set_allele_frequencies(0.2 * np.ones(L), N)

    # run for G generations to measure execution time
    t1=time.time()
    pop.evolve(G)
    t2=time.time()
    
    print "time required for",G,"generations:",round(t2-t1,3),'s'
    exec_time_single_xo.append([L, t2-t1])  # store the execution time
    
exec_time_single_xo=np.array(exec_time_single_xo)


# 3. without recombination
print "\nno recombination"
exec_time_norec = []
for L in range(2,Lmax_single_xo+1):
    print L,"loci, maximum",Lmax_single_xo,

    # set up population
    pop = h.haploid_lowd(L)     # produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = N   # set the population size

    # set and additive fitness function with random coefficients. 
    pop.set_fitness_additive(0.01 * np.random.randn(L))
    
    pop.set_mutation_rates(mu)      # mutation rate

    # initialize the population with N individuals in linkage equilibrium
    pop.set_allele_frequencies(0.2 * np.ones(L), N)

    # evolve for G generation skipping the recombination step
    t1=time.time()
    pop.evolve_norec(G)
    t2=time.time()

    print "time required for",G,"generations:",round(t2-t1,3),'s'
    exec_time_norec.append([L, t2-t1])      # store the execution time

exec_time_norec=np.array(exec_time_norec)


# Plot the execution times of all three cases and
# lines indicating the expected computational complexity
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
