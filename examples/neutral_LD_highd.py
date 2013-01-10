# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       23/08/12
content:    Example of haploid_highd demonstrating the balance between
            recombination and genetic drift in a finite population. The example
            plots the average r^2 (squared correlation coefficient) for loci at
            different distances along the chromosome.
'''
# Import modules (setting the path should not be necessary when the module is
# installed in the PYTHONPATH)
import sys
sys.path.insert(0, '../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h


# specify parameters
N = 500                            # population size
L = 1000                           # number of loci
mu = 0.1 / N                       # mutation rate
r = 10.0 / L / N                   # crossover rate

# set up population
pop = h.haploid_highd(L)             # produce an instance of haploid_lowd with L loci
pop.carrying_capacity = N            # set the steady-state population size
pop.outcrossing_rate=1               # obligate sexual
pop.crossover_rate = r               # crossover rate per locus
pop.mutation_rate=mu                 # mutation rate per locus

# initialize the population with
#    - N/2 individuals with genotypes 00000000...0000
#    - N/2 with the opposite genotype 11111111...1111
pop.set_genotypes([np.zeros(L), np.ones(L)],[N/2, N/2])

# locus pairs for which LD is to be tracked
ld_points = np.arange(5,L-1,100)
locus_pairs = [ [L/2, l1] for l1 in ld_points]

print "Evolve for >> N generations and measure correlations between loci, aka LD."
pop.status()
pop.evolve(10 * N)                      # evolve for 10N generations to equilibrate

pop.evolve(4*N)                     # evolve for 4N to equilibrate

# evolve the population and track linkage disequilibrium (LD)
nsamples = 10000      
LD=[]
rsq = []
for ii in range(nsamples):
    pop.evolve(0.1 * N)             # N / 10 generations between successive samples

    if (ii%100==0):
        print ii, "out of", nsamples, "samples"

    # get allele frequencies
    af = pop.get_allele_frequencies()

    # get LD and r^2
    templd = [pop.get_LD(l1,l2) for l1,l2 in locus_pairs]
    LD.append(templd)
    rsq.append([templd[pi]**2/(af[l1]*(1-af[l1])*af[l2]*(1-af[l2])+1e-10) for pi,(l1,l2) in enumerate(locus_pairs)])
    

# plot the result
plt.figure()
plt.plot(ld_points-L/2, np.mean(rsq,axis=0))
plt.xlabel('Distance on genome')
plt.ylabel(r'$\langle r^2 \rangle$')
plt.title("".join(map(str,[r'$N=',N,r',\,\rho=',r,r',\,\mu=',mu,'$'])))

plt.ion()
plt.show()
