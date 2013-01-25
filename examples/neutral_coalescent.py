import FFPopSim as h
import numpy as np
from matplotlib import pyplot as plt
import random as rd

L=100

pop=h.haploid_highd(L)
pop.outcrossing_rate=0.1
pop.crossover_rate=1.0/pop.L
pop.mutation_rate=5.0/pop.L
pop.carrying_capacity=30

#track the loci 10, 50 and 90
pop.track_locus_genealogy([10,50,90])

#initialize the populations
pop.set_wildtype(pop.carrying_capacity)

pop.status()

#evolve population for several coalescent times
pop.evolve(4*pop.N)

#get tree at locus 10
tree = pop.genealogy.get_tree(10)
print "\nTree statistics:"
print "Branch length:", tree.total_branch_length()
print "External branch length:", tree.external_branch_length()
print "Number of leafs:", len(tree.leafs)
print "Number of internal nodes:", len(tree.nodes)-len(tree.leafs)
print "Time to MRCA:", pop.generation - tree.MRCA.age


#produce a subtree with of a sample of leafs
n = 5
subsample = rd.sample(tree.leafs, n)
sub_tree = tree.create_subtree_from_keys(subsample)
print "\nSubtree with",n,"leafs"
print sub_tree.print_newick()
print "Each tree label is composed of the index if the individual and the size of the clone"

#trees can be exported as a BioPython tree structure
from Bio import Phylo as P
BPtree = tree.to_Biopython_tree()
plt.figure()
P.draw(BPtree)

#in absence of recombination, trees at all three loci are identical. 
#with crossovers, tree decouple with increasing outcrossing rate. 
for locus in pop.genealogy.loci:
	BPtree = pop.genealogy.get_tree(locus).to_Biopython_tree()
	plt.figure()
	plt.title('Tree at locus '+str(locus))
	P.draw(BPtree)

