import FFPopSim as h
import numpy as np
from matplotlib import pyplot as plt
import random as rd

L=1000
Lsel = 200
s = 0.02

pop=h.haploid_highd(L)
pop.outcrossing_rate=0.01
pop.crossover_rate=1.0/pop.L
pop.mutation_rate=0.12/pop.L
pop.carrying_capacity=1e5

fitness = np.zeros(L)
fitness[:Lsel] = np.random.exponential(0.02, size=Lsel)
pop.set_fitness_additive(fitness)

#track the loci 10, 990
pop.track_locus_genealogy([10, 990])

#initialize the populations
pop.set_wildtype(10)

pop.status()

#evolve population, take 3 samples per generation
pop.tree_sample=3
step=5
trajectory = []
for i in range(200//step):
    pop.evolve(step)
    trajectory.append((pop.generation, pop.N, pop.get_fitness_statistics().mean,
                       pop.get_allele_frequencies()[:Lsel].sum(), pop.get_allele_frequencies()[Lsel:].sum()))
    print(trajectory[-1])


tree1 = pop.genealogy.get_tree(10)
tree2 = pop.genealogy.get_tree(990)

sub_tree1 = tree1.create_subtree_from_keys(tree1.sampled_leafs)
sub_tree2 = tree2.create_subtree_from_keys(tree2.sampled_leafs)

from Bio import Phylo as P
plt.figure()
BPtree1 = sub_tree1.to_Biopython_tree()
BPtree1.ladderize()
P.draw(BPtree1, label_func=lambda x:'', axes=plt.subplot(121))

BPtree2 = sub_tree2.to_Biopython_tree()
BPtree2.ladderize()
P.draw(BPtree2, label_func=lambda x:'', axes=plt.subplot(122))
plt.savefig('trees.png')

trajectory = np.array(trajectory)
plt.figure()
plt.plot(trajectory[:,0], trajectory[:,2])
plt.savefig('fitness.png')

plt.figure()
plt.plot(trajectory[:,0], trajectory[:,3])
plt.plot(trajectory[:,0], trajectory[:,4])
plt.savefig('divergence.png')

with open('test.fasta', 'w') as fh:
    fh.write(tree1.print_sequences())

