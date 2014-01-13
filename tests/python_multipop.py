# vim: fdm=indent
'''
author:     Pavel Sagulenko
date:       10/09/13
content:    Test script for the python bindings (multi_population)
'''

# Import module
import sys
sys.path.insert(0, '../pkg/python')
import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as ff
from Bio import Phylo


# Globals
L = 1000    # number of loci
N = 300    # carrying capacity (population size)
locations = 10

# Construct class
pop = ff.multi_population (locations, L)
pop.carrying_capacity = 1000
pop.track_locus_genealogy([3, 60])
pop.set_random_genotype(10)
pop.mutation_rate=1e-3
for i in xrange (100):
    pop.evolve(10)
    print "generation:" + str(i*10)
