import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h
import time

### specify parameters
N=1000                              #Population size
Lmax=13
r = 0.01
mu=0.001
G=1000

exec_time = []
for L in range(2,Lmax):
    t1=time.time()
    ### set up
    pop = h.haploid_lowd(L)             #produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = N           #set the population size
    pop.set_fitness_additive(0.01*np.random.randn(L))     #set and additive fitness function. Note that FFPopSim models fitness landscape 
                                        #in a +/- rather than 0/1 basis, hence the factor 1/2
    pop.set_recombination_rates(r*np.ones(L-1))      #assign the recombination rates
    pop.set_mutation_rates(mu)           #assign the mutation rate
    
    pop.set_allele_frequencies(0.2*np.ones(L))          #initialize the population with N individuals with genotypes 0, that is ----
    pop.evolve(G)                    #run for 10N generations to equilibrate
    
    t2=time.time()
    exec_time.append([L, t2-t1])
    
###repeat without recombination
exec_time_norec = []
for L in range(2,Lmax):
    t1=time.time()
    ### set up
    pop = h.haploid_lowd(L)             #produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = N           #set the population size
    pop.set_fitness_additive(0.01*np.random.randn(L))     #set and additive fitness function. Note that FFPopSim models fitness landscape 
                                        #in a +/- rather than 0/1 basis, hence the factor 1/2
    pop.set_mutation_rates(mu)           #assign the mutation rate
    
    pop.set_allele_frequencies(0.2*np.ones(L))          #initialize the population with N individuals with genotypes 0, that is ----
    pop.evolve_norec(G)                    #run for 10N generations to equilibrate
    
    t2=time.time()
    exec_time_norec.append([L, t2-t1])


exec_time=np.array(exec_time)
exec_time_norec=np.array(exec_time_norec)

plt.figure()
plt.plot(exec_time[:,0], exec_time[:,1],label='with recombination', linestyle='None', marker = 'o')
plt.plot(exec_time[:,0], exec_time[-1,1]/3.0**(Lmax-exec_time[:,0]-1),label=r'$\propto 3^L$')

plt.plot(exec_time_norec[:,0], exec_time_norec[:,1],label='without recombination', linestyle='None', marker = 'x')
plt.plot(exec_time[:,0], exec_time_norec[-1,1]/2.0**(Lmax-exec_time_norec[:,0]-1),label=r'$\propto 2^L$')

plt.plot(exec_time[:,0], exec_time[-1,1]/3.0**(Lmax)*8**(exec_time[:,0]),label=r'$\propto 8^L$')

ax=plt.gca()
ax.set_yscale('log')
plt.xlabel('number of loci')
plt.ylabel('seconds for '+str(G)+' generations')
plt.legend(loc=2)
plt.xlim([1,Lmax])
plt.ylim([0.2*np.min(exec_time_norec[:,1]),10*np.max(exec_time[:,1])])
