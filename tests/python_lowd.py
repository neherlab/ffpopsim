# vim: fdm=indent
'''
author:     Fabio Zanini
date:       14/05/12
content:    Test script for the python bindings to the low-dimensional
            simulation
'''

# Import module
import sys
sys.path.append('../pkg/python')
import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h

# Construct class
pop = h.haploid_lowd(4)

# Test initialization
#pop.set_allele_frequencies([0,0.3,0.6,0.9], 1000)
pop.set_genotypes([1,2],[400,800])

# Test setting the recombination/mutation rates
pop.set_recombination_rates([0.01, 0.03, 0.02], h.SINGLE_CROSSOVER)
pop.set_mutation_rates([0.003,0.002,0.004,0.005],
                       [0.006,0.004,0.008,0.010])

# Test getting the mutation rate
print pop.get_mutation_rates(direction=0)
print pop.get_mutation_rates(direction=1)

# Test setting / getting fitness
pop.set_fitness_additive([0.02,0.03,0.04,0.02])
pop.get_fitnesses()

# Test allele frequency readout
print pop.get_allele_frequencies()

# Test evolution
gens = 100
from time import time as ti
t0 = ti()
pop.evolve(gens)
t1 = ti()
print 'Time for evolving the population for '+str(gens)+' generations: {:1.1f} s'.format(t1-t0)

# Print population size
print pop.N

# Test divergence / diversity statistics
print pop.get_divergence_statistics()
print pop.get_diversity_statistics()

# Plot histograms
plt.ion()
pop.plot_fitness_histogram()
pop.plot_divergence_histogram(color='r')
pop.plot_diversity_histogram(color='g')
