# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       24/08/12
content:    Example of haploid_lowd on the algorithm runtime complexity.
            The same series of simulations are performed:
                - for a human-like population,
                - for an HIV-like (viral) population,
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


# 1. specify parameters for a human population, per MB of DNA
r = 1e-8                       # crossover rate
mu = 1e-8                      # mutation rate
G = 100                        # generations
NS = 1000
exec_time_MB = {}
Llist_MB = [1e5, 3e5, 1e6]
Nlist = [100, 300, 1000, 3000, 10000, 30000]

for L in Llist_MB:
    exec_time = []
    for N in Nlist:
        # set up population
        pop = h.haploid_highd(L)    # produce an instance of haploid_lowd with L loci
        pop.carrying_capacity = N   # set the population size

        pop.recombination_model = h.CROSSOVERS  # set recombination model
        pop.outcrossing_rate = 1.0              # obligate sexual
        pop.crossover_rate = r                  # crossover rates
        pop.mutation_rate = mu                  # mutation rate (per locus)
        
        pop.set_wildtype(N)                     # set a wildtype population of size N

        pop.evolve(1.0 / (L * (mu + r)))        # evolve until equilibrium
        
        # run for G generations to measure execution time    
        t1=time.time()
        pop.evolve(G)
        t2=time.time()
        print L, "bases, human like, population size:", N,
        print "time required for", G, "generations:",
        print round(t2 - t1, 3), "seconds"
        
        # set additive fitness landscape and repeat
        selection_coefficients = np.zeros(L)
        for locus in xrange(0, int(L), int(L/NS)):
            selection_coefficients[locus] = 0.01 * ((np.random.rand()-0.98)>0)
        pop.set_fitness_additive(selection_coefficients)

        # run for G generations to measure execution time (with fitness)    
        t3=time.time()
        pop.evolve(G)
        t4=time.time()

        print L, "bases, human like, population size:", N,
        print "time required for", G, "generations with selection:",
        print round(t4 - t3, 3), "seconds"

        exec_time.append([N, t2 - t1, t4 - t3])    # store the execution time
    
    exec_time_MB[L] = np.array(exec_time)

# plot human curves
cols=['r','g','b','c','m','k']
plt.figure()
for ii,L in enumerate(Llist_MB):
    plt.plot(exec_time_MB[L][:,0],
             exec_time_MB[L][:,1],
             label = r'human like, $L=10^{'+str(round(np.log10(L),2))+'}$',
             c=cols[ii])

    plt.plot(exec_time_MB[L][:,0],
             exec_time_MB[L][:,2],
             ls='--',
             c=cols[ii])


# 2. specify parameters for a virus-like genome
r = 1e-3            # recombination rate
coinf = 0.01        # coinfection rate  == outcrossing rate  
mu = 1e-5           # mutation rate
G = 100             # generations

exec_time_virus = {}
Nlist = [1000, 3000, 10000, 30000, 100000, 300000, 1000000]
Llist_virus = [1e3, 3e3, 1e4]
for L in Llist_virus:
    exec_time = []
    for N in Nlist:

        # set up population
        pop = h.haploid_highd(L)        # produce an instance of haploid_lowd with L loci
        pop.carrying_capacity = N       # set the population size
        pop.crossover_rate = r          # crossover rates
        pop.mutation_rate = mu          # mutation rate (per locus)
        pop.outcrossing_rate = coinf    # outcrossing rate
        
        pop.set_wildtype(N)             # set a wildtype population of size N

        # run for G generations to measure execution time
        t1=time.time()
        pop.evolve(G)
        t2=time.time()
        print "Population size:", N, "genome length:",L
        print "time required for", G, "generations:",
        print round(t2 - t1, 3), "seconds"

        # set additive fitness landscape and repeat
        selection_coefficients = np.zeros(L)
        for locus in xrange(0, int(L), int(L/NS)):
            selection_coefficients[locus] = 0.01 * ((np.random.rand()-0.98)>0)
        pop.set_fitness_additive(selection_coefficients)

        # run for G generations to measure execution time (with fitness)    
        t3=time.time()
        pop.evolve(G)    
        t4=time.time()
        print "Population size:",N,"time required for", G ,"generations with selection:",round(t4-t3,3),"seconds"

        exec_time.append([N, t2-t1,t4-t3])    # store the execution time
        
    exec_time_virus[L]=np.array(exec_time)


# plot viral curves
for ii,L in enumerate(Llist_virus):
    plt.plot(exec_time_virus[L][:,0],
             exec_time_virus[L][:,1],
             label =r'virus like, $L=10^{'+str(round(np.log10(L),2))+'}$',
             c=cols[ii+3])

    plt.plot(exec_time_virus[L][:,0],
             exec_time_virus[L][:,2],
             ls='--',
             c=cols[ii+3])

ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.legend(loc=4)
plt.xlabel('Population size')
plt.ylabel('seconds for '+str(G)+' generations')

plt.ion()
plt.show()
