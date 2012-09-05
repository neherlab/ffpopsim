'''
author:     Richard Neher
date:       11/07/12
content:    Example on genetic drift using haploid_highd
'''
# Import module (setting the path should not be necessary when the module is installed in the python path
import sys
sys.path.insert(0,'../pkg/python')

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import FFPopSim as h

#simulate 256 loci
L=256

### set up
pop = h.haploid_highd(L)                        #produce an instance of haploid_lowd with L loci
pop.carrying_capacity = 5000                    #set the average population size to 50000
pop.outcrossing_rate = 1                        #make the species obligate outcrossing
pop.crossover_rate = 0.02/pop.L                 #set the crossover rate of the segment to 2 centimorgans
pop.mutation_rate = 0.1/pop.carrying_capacity   #per locus mutation rate equal to 0.1/N


initial_allele_frequencies = 0.5*np.ones(pop.L) #define some initial allele frequencies

#initialize the population in LD with the specified allele frequencies
pop.set_allele_frequencies(initial_allele_frequencies, pop.carrying_capacity)

#make lists to store allele frequencies and time points
allele_frequencies = [pop.get_allele_frequencies()] 
tp = [pop.generation]

#evolve for 2000 generations and track the allele frequencies
maxgen = 2000
while pop.generation<maxgen:
    pop.evolve(10)                                          #procede 10 generations
    allele_frequencies.append(pop.get_allele_frequencies()) #save the allele frequencies
    tp.append(pop.generation)                               #and the associated generation

#convert to an array to enable slicing
allele_frequencies=np.array(allele_frequencies)

#plot the result
plt.figure()
for locus in xrange(5,pop.L,50):        #plot a few neutral trajectories
    plt.plot(tp, allele_frequencies[:,locus], c=cm.cool(locus), lw=2)

plt.title('Genetic Drift')
plt.xlabel('Time [generations]')
plt.ylabel('Allele frequencies')


