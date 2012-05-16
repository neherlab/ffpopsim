# vim: fdm=indent
'''
author:     Fabio Zanini
date:       14/05/12
content:    Test script for the python bindings to the low-dimensional
            simulation
'''

# Import module
import numpy as np
import matplotlib.pyplot as plt
import PopGenLib as h

# Construct class
c = h.haploid_gt_dis()
c.set_up(4, 1000)

## Test initialization
#c.init_frequencies([0,0.3,0.6,0.9])
#c.init_genotypes([1,2],[0.4,0.8])
#
## Test setting the recombination/mutation rates
#c.set_recombination_rates([0.01, 0.03, 0.02])
#c.set_mutation_rate([[0.003,0.002,0.004,0.005],
#                     [0.006,0.004,0.008,0.010]])

# Test getting the mutation rate
print c.get_mutation_rate(direction=1)

## Test setting / getting fitness
#c.set_fitness_additive([0.02,0.03,0.04,0.02])
#c.get_fitnesses()
#
## Test allele frequency readout
#print c.get_allele_frequencies()
#
## Test evolution
#gens = 100
#from time import time as ti
#t0 = ti()
#c.evolve(gens)
#t1 = ti()
#print 'Time for evolving the population for '+str(gens)+' generations: {:1.1f} s'.format(t1-t0)
#
## Test divergence / diversity statistics
#print c.get_divergence_statistics()
#print c.get_diversity_statistics()
#
### Plot histograms
#plt.ion()
#c.plot_fitness_histogram()
#c.plot_divergence_histogram(color='r')
#c.plot_diversity_histogram(color='g')
