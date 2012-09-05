# vim: fdm=indent
'''
author:     Richard Neher
date:       23/08/12
content:    Example of haploid_highd demonstrating the balance between recombination and 
            genetic drift in a finite population. The example plots the average r^2 for loci
            at different distances along the chromosome.
'''
# Import module (setting the path should not be necessary when the module is installed in the python path
import sys
sys.path.insert(0,'../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h


# specify parameters
N = 500                            # Population size
L = 1000                           # number of loci
mu = 0.1/N                            # mutation rate
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

### equilibrate
pop.evolve(4*N)

#### sample allele frequencies and store in array allele_frequencies
nsamples = 10000                    #increase for better statistics      

### propagate the population and periodically calculate LD 
LD=[]
rsq = []
for ii in range(nsamples):
    pop.evolve(0.1*N)               #N/10 generations between successive samples
    if (ii%100==0): print ii, "out of", nsamples 
    af=pop.get_allele_frequencies()
    templd = [pop.get_LD(l1,l2) for l1,l2 in locus_pairs]
    LD.append(templd)
    rsq.append([templd[pi]**2/(af[l1]*(1-af[l1])*af[l2]*(1-af[l2])+1e-10) for pi,(l1,l2) in enumerate(locus_pairs)])
    

#### plot the result
plt.figure()
plt.plot(ld_points-L/2, np.mean(rsq,axis=0))
plt.xlabel('Distance on genome')
plt.ylabel(r'$\langle r^2 \rangle$')
plt.title("".join(map(str,[r'$N=',N,r',\,\rho=',r,r',\,\mu=',mu,'$'])))
