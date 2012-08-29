# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/04/12
content:    Test script for the python bindings
'''

# Import module
import sys
sys.path.append('../pkg/python')
import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h

# Construct class
pop = h.hivpopulation(1000)

# Test I/O fitness landscapes
pop.set_replication_landscape(lethal_fraction=0.05,
                            number_valleys=0)
pop.read_replication_coefficients('hiv_model.dat')
rep = pop.get_replication_additive()
rep[np.random.random(10000) > 0.5] = -0.1
pop.set_replication_additive(rep)

# Show the additive part of the fitness landscape
print pop.get_trait_additive()

# Test population initialization
pop.set_allele_frequencies([0.3] * h.HIVGENOME, 1000)

# Test allele frequency readout
print np.max(pop.get_allele_frequency(4))

# Test evolution
from time import time as ti
t0 = ti()
pop.evolve(30)
t1 = ti()
print 'Time for evolving HIV for 30 generations: {:1.1f} s'.format(t1-t0)

# Write genotypes
pop.write_genotypes('test.txt', 100)
pop.write_genotypes_compressed('test.npz', 100)

# Plot histograms
plt.ion()
pop.plot_fitness_histogram()
pop.plot_divergence_histogram(color='r')
pop.plot_diversity_histogram(color='g')

# Test treatment changes
pop.treatment = 0.4
print pop.treatment
