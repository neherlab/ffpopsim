# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/08/12
content:    Test script for the python bindings (haploid_highd)
'''

# Import module
import sys
sys.path.insert(0, '../pkg/python')
import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h
from Bio import Phylo


# Globals
L = 1000    # number of loci
N = 300    # population size


# Construct class
pop = h.haploid_highd(L, all_polymorphic=False)

# Set a growth rate explicitely
pop.growth_rate = 2.5

# Start tracking genealogy of two loci
pop.track_locus_genealogy([3, 60])

# Test fitness landscapes
rep = np.zeros(L)
rep[np.random.random(L) > 0.5] = -0.1
pop.set_trait_additive(rep)

# Show the additive part of the fitness landscape
print pop.get_trait_additive()

# Test population initialization
pop.track_locus_genealogy([3,6])
pop.set_wildtype(N)
#pop.set_allele_frequencies([0.3] * L, N)
pop.mutation_rate = 1e-5
pop.outcrossing_rate = 1e-2
pop.crossover_rate = 1e-3

# Test allele frequency readout
print np.max(pop.get_allele_frequency(4))

# Test evolution
from time import time as ti
t0 = ti()
pop.evolve(30)
t1 = ti()
print 'Time for evolving population for 30 generations: {:1.1f} s'.format(t1-t0)

## Write genotypes
#pop.write_genotypes('test.txt', 100)
#pop.write_genotypes_compressed('test.npz', 100)

## Plot histograms
#plt.ion()
#pop.plot_fitness_histogram()
#pop.plot_divergence_histogram(color='r')
#pop.plot_diversity_histogram(color='g')

# Look at the genealogy
print pop.genealogy
tree = pop.genealogy.get_tree(3)
subtree = tree.create_subtree_from_keys([tree.leafs[1], tree.leafs[2]])

# Convert the tree into Biopython format and plot it
treeBio = tree.to_Biopython_tree()
Phylo.draw(treeBio)
