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
N = 500                            # Population size
L = 1000                           # number of loci
mu = 0.5/N                            # mutation rate
r = 10.0/L/N                           # crossover rate

### set up
pop = h.haploid_highd(L)             # produce an instance of haploid_lowd with L loci
pop.carrying_capacity = N            # set the steady-state population size
pop.outcrossing_rate=1               # assign the outcrossing rate (equals 1 be default)
pop.crossover_rate = r               # assign the crossover rate per locus
pop.mutation_rate=mu                 # assign the mutation rate per locus

# initialize the population with N/2 individuals with genotypes 00000000...0000
# and N/2 with the opposite genotype 111111111111....111
pop.set_genotypes([np.zeros(L), np.ones(L)],[N/2, N/2])

#locus pairs for which LD is to be tracked
ld_points = np.arange(5,L-1,100)
locus_pairs = [ [L/2, l1] for l1 in ld_points]

pop.evolve(4*N)

#### sample allele frequencies and store in array allele_frequencies
nsamples = 10000                    #increase for better statistics      

LD=[]
rsq = []
for ii in range(nsamples):
    pop.evolve(0.1*N)               #N/10 generations between successive samples
    if (ii%100==0): print ii, "out of", nsamples 
    af=pop.get_allele_frequencies()
    templd = [pop.get_LD(l1,l2) for l1,l2 in locus_pairs]
    LD.append(templd)
    rsq.append([templd[pi]**2/(af[l1]*(1-af[l1])*af[l2]*(1-af[l2])+1e-10) for pi,(l1,l2) in enumerate(locus_pairs)])
    
LD=np.array(LD);
rsq=np.array(rsq);



plt.figure()
plt.plot(ld_points, np.mean(rsq,axis=0))
rho = 2*r*N*abs(L/2-ld_points)
plt.plot(ld_points, (rho+10)/(rho**2+13*rho+22))
