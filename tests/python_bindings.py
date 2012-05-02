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

# Test population initialization
# FIXME: tends to freeze (fittest clone explodes!)
print c.init_genotypes(np.zeros(h.HIVGENOME) + 0.3, 1000)

## Plot fitness distribution
#ncl = c.get_number_of_clones()
#fit = c.get_fitnesses(ncl);
#plt.ion()
#plt.figure()
#plt.hist(fit)

# Test allele frequency readout
print np.max(c.get_allele_frequency(4))

# Test evolution
#c.evolve(1)

## Plot fitness distribution
#ncl = c.get_number_of_clones()
#fit = c.get_fitnesses(ncl);
#plt.figure()
#plt.hist(fit)

# Test treatment changes
c.set_treatment(0.4)
c.get_treatment()


