# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/04/12
content:    Test script for the python bindings
'''

# Import module
import numpy as np
import matplotlib.pyplot as plt
import PopGenLib as h

# Construct class
c = h.hivpopulation()
c.set_up(1000)

# Test I/O fitness landscapes
c.read_selection_coefficients('hiv_model.dat')

## Test population initialization
#print c.init_frequencies(np.zeros(h.HIVGENOME) + 0.3, 1000)

# Test allele frequency readout
print np.max(c.get_allele_frequency(4))

# Test evolution
from time import time as ti
t0 = ti()
c.evolve(100)
t1 = ti()
print 'Time for evolving HIV for 100 generations: {:1.1f} s'.format(t1-t0)

# Write genotypes
c.write_genotypes('test.txt', 100)
c.write_genotypes_compressed('test.npz', 100)

# Plot histograms
plt.ion()
c.plot_fitness_histogram()
c.plot_divergence_histogram(color='r')
c.plot_diversity_histogram(color='g')

# Test treatment changes
c.set_treatment(0.4)
c.get_treatment()
