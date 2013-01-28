import FFPopSim as h
import numpy as np
from matplotlib import pyplot as plt
import random as rd
from Bio import Phylo

print "This script is meant to illustrate and explore the effect of\n\
positive selection on genealogies in asexual and sexual populations. \n\n\
Simulations are performed using an infinite sites model with L segregating\n\
sites at which mutations with identical beneficial effect are injected.\n\n"

#suggested values
#neutral asexual:	N=100 	s=0.00001	r=0.0
#selected asexual: 	N=10000	s=0.01		r=0.0
#selected sexual: 	N=1000 	s=0.01		r=1.0

L = 1000   	#number of segregating sites
s = 1e-2 	#single site effect
N = 10000 	#population size
r = 0.0  	#outcrossing rate

sample_size=30	#number of individuals whose genealogy is looked at
nsamples = 3	#number of trees
burnin = 2000 	#either ~5*N or 5/s, depending on whether coalescence is dominated by drift or draft
dt = 1000 		#time between samples

#set up population, switch on infinite sites mode
pop=h.haploid_highd(L, all_polymorphic=True)

#set the population size via the carrying capacity
pop.carrying_capacity= N

#set the crossover rate, outcrossing_rate and recombination model
pop.outcrossing_rate = r
pop.recombination_model = h.CROSSOVERS
pop.crossover_rate = 1.0/pop.L

#set the effect sizes of the mutations that are injected (the same at each site in this case)
pop.set_fitness_additive(np.ones(L)*s)

#track the genealogy at a central locus L/2 (which one doesn't matter in the asexual case)
pop.track_locus_genealogy([L/2])

#initialize the populations
pop.set_wildtype(pop.carrying_capacity)

print "Population parameters:"
pop.status()

#burn in
print "\nEquilibrate:"
while pop.generation<burnin:
	print "Burn in: at", pop.generation, "out of", burnin, "generations"
	pop.evolve(100)


print "\nPlot coalescent trees:"
fig=plt.figure(figsize=(7,10))
fig.suptitle("".join(map(str,['N=',N,'  r=',r,'  L=',L, '  s=',s])), fontsize=18)
for si in xrange(nsamples):
	print "sample",si,"out of",nsamples
	#evolve a while before sampling the next tree
	pop.evolve(dt)

	#draw a sample from the population, convert its genealogy to a BioPython tree object and plot
	tree = pop.genealogy.get_tree(L/2)
	subtree = tree.create_subtree_from_keys(rd.sample(tree.leafs,sample_size)).to_Biopython_tree()
	subtree.ladderize()
	plt.subplot(3,1,si+1)
	Phylo.draw(subtree,label_func=lambda x:"")
	plt.draw()

plt.savefig("".join(map(str,['tree_', 'N=',N,'_r=',r,'_L=',L, '_s=',s,'.pdf'])))

