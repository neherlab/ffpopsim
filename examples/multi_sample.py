import FFPopSim as h
import numpy as np
from matplotlib import pyplot as plt
import random as rd

L=100

pop=h.haploid_highd(L)
pop.outcrossing_rate=0.0
pop.crossover_rate=1.0/pop.L
pop.mutation_rate=0.01/pop.L
pop.carrying_capacity=50

#track the loci 10, 50 and 90
pop.track_locus_genealogy([10,50,90])

#initialize the populations
pop.set_wildtype(pop.carrying_capacity)

pop.status()

#evolve population for several coalescent times to burn in
pop.evolve(4*pop.N)

for sample in range(50):
    # set sample size to 10, evolve on generation to take a sample
    pop.tree_sample=5
    pop.evolve(1)
    # set sample to 0, evolve for N/3
    pop.tree_sample=0
    pop.evolve(pop.N//2)

pop.tree_sample=5
pop.evolve(1)

tree = pop.genealogy.get_tree(10)

sub_tree = tree.create_subtree_from_keys(tree.sampled_leafs)

from Bio import Phylo as P
BPtree = sub_tree.to_Biopython_tree()
BPtree.ladderize()
plt.figure()
P.draw(BPtree, label_func=lambda x:'')

plt.savefig('tree.png')

with open('test.fasta', 'w') as fh:
    fh.write(tree.print_sequences())

