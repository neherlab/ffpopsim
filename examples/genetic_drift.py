'''
author:     Richard Neher, Fabio Zanini
date:       11/07/12
content:    Example on genetic drift using haploid_highd
'''
# Import modules (setting the path should not be necessary when the module is
# installed in the PYTHONPATH)
import sys
sys.path.insert(0, '../pkg/python')

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import FFPopSim as h


# specify parameters
L = 256                                         # simulate 256 loci

# set up population
pop = h.haploid_highd(L)                        # produce an instance of haploid_highd with L loci
pop.carrying_capacity = 50000                   # set the average population size to 50000
pop.outcrossing_rate = 1                        # make the species obligate outcrossing
pop.crossover_rate = 0.02 / pop.L               # set the crossover rate of the segment to 2 centimorgans
pop.mutation_rate = 0.1 / pop.carrying_capacity # per locus mutation rate equal to 0.1/N


# initialize the population in linkage equilibrium with the specified allele frequencies
initial_allele_frequencies = 0.5*np.ones(pop.L) # define some initial allele frequencies as 1/2
pop.set_allele_frequencies(initial_allele_frequencies, pop.carrying_capacity)

# evolve for 2000 generations and track the allele frequencies
maxgen = 2000
allele_frequencies = [pop.get_allele_frequencies()] 
tp = [pop.generation]

print "Illustrate genetic drift on allele frequency trajectories."
pop.status()    #print status message
while pop.generation < maxgen:
    if (pop.generation%(maxgen/10)==0): print pop.generation,"out of",maxgen, "generations"
    pop.evolve(10)

    # save allele frequencies and time
    allele_frequencies.append(pop.get_allele_frequencies())
    tp.append(pop.generation)

# convert to an array to enable slicing
allele_frequencies = np.array(allele_frequencies)

# plot the result
plt.figure()
for locus in xrange(5,pop.L,50):                # plot a few neutral trajectories
    plt.plot(tp, allele_frequencies[:,locus], c=cm.cool(locus), lw=2)

plt.title('Genetic Drift')
plt.xlabel('Time [generations]')
plt.ylabel('Allele frequencies')

plt.ion()
plt.show()

