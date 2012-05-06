#!/usr/bin/env python2
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/04/12
content:    Test script for the python bindings
'''

# Import module
import numpy as np
import matplotlib.pyplot as plt
import hivpython as h

# Construct class
c = h.hivpython()
c.set_up(1000)

# Test I/O fitness landscapes
c.read_selection_coefficients('hiv_model.dat')

## Test population initialization
#print c.init_genotypes(np.zeros(h.HIVGENOME) + 0.3, 1000)

# Test allele frequency readout
print np.max(c.get_allele_frequency(4))

# Test evolution
c.evolve(10)

# Write genotypes
#c.write_genotypes('test.txt', 100)
c.write_genotypes_compressed('test.npz', 100)

## Plot histograms
#plt.ion()
#c.plot_fitness_histogram()
#c.plot_divergence_histogram(color='r')
#
## Test treatment changes
#c.set_treatment(0.4)
#c.get_treatment()


