import FFPopSim as h
import numpy as np
from matplotlib import pyplot as plt
import random as rd

def get_mut_count_subtree(R):
	if R.is_terminal():
		return 1,[]
	else:
		tmp_mut_count = []
		sub_tree_size = 0
		for C in R.clades:
			csize,cmut  = get_mut_count_subtree(C)
			sub_tree_size += csize
			tmp_mut_count.extend(cmut)
			tmp_mut_count.append([csize,C.branch_length])
		return sub_tree_size, tmp_mut_count


def get_SFS(T):
	sample_size,SFS = get_mut_count_subtree(T.root)
	SFS = np.asarray(SFS)
	SFS[:,0]/=sample_size
	return sample_size,SFS

L=100

pop=h.haploid_highd(L)
pop.outcrossing_rate=1.0
pop.crossover_rate=1.0/pop.L
pop.mutation_rate=0.0/pop.L
pop.carrying_capacity=100

#track the loci 10, 50 and 90
pop.track_locus_genealogy([10,50,90])

#initialize the populations
pop.set_wildtype(pop.carrying_capacity)

pop.status()

#evolve population for several coalescent times
pop.evolve(4*pop.N)

nsamples = 100
nloci = len(pop.genealogy.loci)
SFS = []
for si in xrange(nsamples):
	pop.evolve(pop.N)
	print "sample", si, "out of", nsamples
	for locus in pop.genealogy.loci:
		BPtree = pop.genealogy.get_tree(locus).to_Biopython_tree()
		sample_size,tmpSFS = get_SFS(BPtree)
		SFS.extend(tmpSFS.tolist())
	
SFS=np.asarray(SFS)
y,x = np.histogram(SFS[:,0], weights = SFS[:,1], bins=np.logspace(-1.5,0,11))
bincenters = 0.5*(x[1:]+x[:-1])
dx= x[1:]-x[:-1]
plt.plot(bincenters, y/nsamples/nloci/dx, label='simulation')
plt.plot(bincenters, 2*pop.N/bincenters, label = 'coalescent')
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.ylabel('SFS')
plt.xlabel('derived allele frequency')
plt.legend()
