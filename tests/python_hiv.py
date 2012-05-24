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
c = h.hivpopulation(1000)

# Test I/O fitness landscapes
c.set_replication_landscape(lethal_fraction=0.05,
                            number_valleys=3)
#c.read_replication_coefficients('hiv_model.dat')

# Show the additive part of the fitness landscape
print c.get_additive_trait()

# Test population initialization
c.set_allele_frequencies([0.3] * h.HIVGENOME, 1000)

# Test allele frequency readout
print np.max(c.get_allele_frequency(4))

# Test evolution
from time import time as ti
t0 = ti()
c.evolve(30)
t1 = ti()
print 'Time for evolving HIV for 30 generations: {:1.1f} s'.format(t1-t0)

# Write genotypes
c.write_genotypes('test.txt', 100)
c.write_genotypes_compressed('test.npz', 100)

# Plot histograms
plt.ion()
c.plot_fitness_histogram()
c.plot_divergence_histogram(color='r')
c.plot_diversity_histogram(color='g')

# Test treatment changes
c.treatment = 0.4
print c.treatment
