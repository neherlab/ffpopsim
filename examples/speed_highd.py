# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       24/08/12
content:    Example of haploid_lowd on the algorithm complexity
'''
# Import module
import sys
sys.path.insert(0,'../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h
import time

# specify parameters for a MB of DNA
r = 1e-8                        # crossover rate
mu = 1e-8                      # Mutation rate
G = 100                         # Generations

exec_time_MB = {}
Llist_MB = [1e5, 3e5,1e6]
Nlist = [100,300,1000,3000,10000,30000]
for L in Llist_MB:
    exec_time = []
    for N in Nlist:
        ### set up
        pop = h.haploid_highd(L)    # produce an instance of haploid_lowd with L loci
        pop.carrying_capacity = N   # set the population size
        pop.crossover_rate = r      # assign the crossover rates
        pop.mutation_rate= mu       # assign the mutation rate (per locus)
        pop.recombination_model = h.CROSSOVERS
        pop.outcrossing_rate = 1.0  # make outcrossing obligatory
        
        pop.set_wildtype(N)
        pop.evolve(1.0/(L*(mu+r)))
        
        t1=time.time()
        pop.evolve(G)                   # run for G generations to measure execution time    
        t2=time.time()
        print L, "bases, human like, population size:",N,"time required for", G ,"generations:",round(t2-t1,3),"seconds"
        exec_time.append([N, t2-t1])    # store the execution time
    
    exec_time_MB[L]=np.array(exec_time)


# specify parameters for a virus like genome
r = 1e-3            # Recombination rate
coinf = 0.01        # coinfection rate  == outcrossing rate  
mu = 1e-5           # Mutation rate
G = 100             # Generations

exec_time_virus = {}
Nlist = [1000,3000,10000,30000,100000,300000,1000000]
Llist_virus = [1e3, 3e3,1e4]
for L in Llist_virus:
    exec_time = []
    for N in Nlist:
        ### set up
        pop = h.haploid_highd(L)    # produce an instance of haploid_lowd with L loci
        pop.carrying_capacity = N   # set the population size
        pop.crossover_rate = r      # assign the crossover rates
        pop.mutation_rate= mu       # assign the mutation rate (per locus)
        pop.outcrossing_rate = coinf  # make outcrossing obligatory
        
        pop.set_wildtype(N)
        t1=time.time()
        pop.evolve(G)                   # run for G generations to measure execution time    
    
        t2=time.time()
        print "Population size:",N,"time required for", G ,"generations:",round(t2-t1,3),"seconds"
        exec_time.append([N, t2-t1])    # store the execution time
        
    exec_time_virus[L]=np.array(exec_time)


plt.figure()
for L in Llist_MB:
    plt.plot(exec_time_MB[L][:,0], exec_time_MB[L][:,1], label = 'Human like, L='+str(L))

for L in Llist_virus:
    plt.plot(exec_time_virus[L][:,0], exec_time_virus[L][:,1], label = 'virus like, L='+str(L))

ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.legend(loc=2)
plt.xlabel('Population size')
plt.ylabel('seconds for '+str(G)+' generations')
