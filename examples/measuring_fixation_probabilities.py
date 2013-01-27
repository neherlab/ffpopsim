import FFPopSim as h
import numpy as np
from matplotlib import pyplot as plt
import random as rd

print 'FFPopSim uses a finite sites implementation with a genome \n\
of fixed length L. One can, however, specify the number of segregating sites \n\
and FFPopSim will keep these sites polymorphic by introducing a mutation \n\
in a random individual as soon as the previous polymorphism at this site \n\
disappears or fixes.'

L=1000

pop=h.haploid_highd(L)
pop.outcrossing_rate=0.1
pop.crossover_rate=1.0/pop.L
#switch on infinite sites mode
pop.all_polymorphic=True
pop.mutation_rate=0
#mutation rate needs to be zero in this case. 
pop.carrying_capacity=1000
pop.set_fitness_additive(np.ones(L)*0.01)
#initialize the populations
pop.set_wildtype(pop.carrying_capacity)

pop.status()

#evolve population for several coalescent times
nsamples=100
#bins = np.logspace(-2.5,0,21)
bins = np.exp(np.linspace(-3*np.log(10), 3*np.log(10),21))
bins = bins/(1+bins)
bins[0]=0; bins[-1]=1
bc= 0.5*(bins[:-1]+bins[1:])
dx = bins[1:]-bins[:-1]
SFS=np.zeros(len(bc))
for si in xrange(nsamples):
	print "sample",si,"out of",nsamples
	pop.evolve(100)
	af=pop.get_derived_allele_frequencies()
	y,x = np.histogram(af,bins=bins)
	SFS+=y

plt.plot(np.log(bc/(1-bc)), SFS/dx)
ax=plt.gca()
ax.set_yscale('log')
#ax.set_xscale('log')
plt.xlabel('derived allele frequency')
plt.ylabel('site frequency spectrum')
plt.show()





