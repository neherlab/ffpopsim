import FFPopSim as h
import numpy as np
from matplotlib import pyplot as plt
import random as rd

print 'FFPopSim uses a finite sites implementation with a genome \n\
of fixed length L. One can, however, specify the number of segregating sites \n\
and FFPopSim will keep these sites polymorphic by introducing a mutation \n\
in a random individual as soon as the previous polymorphism at this site \n\
disappears or fixes.'

L=100

pop=h.haploid_highd(L)
pop.outcrossing_rate=0.1
pop.crossover_rate=1.0/pop.L
#switch on infinite sites mode
pop.all_polymorphic=True
pop.mutation_rate=0
#mutation rate needs to be zero in this case. 
pop.carrying_capacity=100

#track the loci 10, 50 and 90
pop.track_locus_genealogy([10,50,90])

#initialize the populations
pop.set_wildtype(pop.carrying_capacity)

pop.status()

#evolve population for several coalescent times
nsamples=100
bins = np.logspace(-1.5,0,11)
bc= 0.5*(bins[:-1]+bins[1:])
dx = bins[1:]-bins[:-1]
SFS=np.zeros(len(bc))
for si in xrange(nsamples):
	pop.evolve(100)
	af=pop.get_derived_allele_frequencies()
	y,x = np.histogram(af,bins=bins)
	SFS+=y

plt.plot(bc, SFS/dx)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('derived allele frequency')
plt.ylabel('site frequency spectrum')
plt.show()





