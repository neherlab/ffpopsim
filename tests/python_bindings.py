#!/usr/bin/env python2
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/04/12
content:    Test script for the python bindings
'''

# Import module
import numpy as np
import hivpython as h

# Construct class
c = h.hivpython()
c.set_up(1000)

# Test I/O fitness landscapes
c.read_selection_coefficients('hiv_model.dat')

# Test evolution
c.evolve(20)

# Test treatment changes
c.set_treatment(0.4)
c.get_treatment()

# Test allele frequency readout
print np.max(c.get_allele_frequency(4))

