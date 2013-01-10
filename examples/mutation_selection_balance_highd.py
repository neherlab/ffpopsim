'''
author:     Richard Neher, Fabio Zanini
date:       11/07/12
content:    Example on the steady state distribution of allele frequency in a 
            balance between mutation and genetic drift using haploid_highd.
'''
# Import modules (setting the path should not be necessary when the module is
# installed in the PYTHONPATH)
import sys
sys.path.insert(0, '../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h


# specify parameters
N = 500                             # population size
L = 64                              # number of loci
s = np.linspace(-2 ,2, L) / N       # additive selection coefficients for L loci, scaled to N
mu = 0.5 / N                        # mutation rate, scaled to N
r  = 50.0 / N / L                   # recombination rate for each interval between loci

# set up population
pop = h.haploid_highd(L)            # produce an instance of haploid_highd with L loci

# set and additive fitness function. Note that FFPopSim models fitness landscape
# in a +/- rather than 0/1 basis, hence the factor 1/2
pop.set_fitness_additive(0.5 * s)

pop.mutation_rate = mu                   # mutation rate
pop.recombination_model = h.CROSSOVERS   # recombination model
pop.outcrossing_rate = 1                 # obligate sexual
pop.crossover_rate = r                   # crossover rate

# initialize population in linkage equilibrium with frequencies 0.5 and size N
pop.carrying_capacity = N           # set the population size
pop.set_allele_frequencies(0.5 * np.ones(L), N)

print "Evolve for >> N generations and compare allele frequency distributions \nto expectations from diffusion theory."
pop.status()
pop.evolve(10 * N)                      # evolve for 10N generations to equilibrate

# evolve and sample allele_frequencies
nsamples = 10000
allele_frequencies = np.zeros((nsamples,L))
for ii in range(nsamples):
    pop.evolve(0.1 * N)                 # N / 10 generations between successive samples

    # print output every 100 generations
    if (ii % 100 == 0):
        print ii, "out of", nsamples, ". Population size: ", pop.population_size, "Number of clones", pop.number_of_clones

    # get allele frequencies
    allele_frequencies[ii,:] = pop.get_allele_frequencies()


# prepare allele frequency histogram
af_bins = np.linspace(0,1,26)           # bins for histogram
af_bins[0] -= 1.0/N
af_bins[-1] += 1.0/N
bc = 0.5*(af_bins[1:]+af_bins[:-1])     # bin centers for plotting

# plot results
plt.figure()
for locus in range(L):

    # make histogram
    y,x = np.histogram(allele_frequencies[:,locus], bins=af_bins, normed='True')

    # plot
    plt.plot(bc, y, color=plt.cm.jet(locus*4))

plt.title('Comparison to diffusion theory for $rN='+str(r*N)+'$, $\mu N='+str(mu*N)+'$, $N='+str(N)+'$')
plt.text(0.3,3,"Color indicates selection coefficient \nfrom Ns=-2..2")
plt.xlabel(r'Allele frequency $\nu$')
plt.ylabel(r'Allele frequency distribution $f(\nu)$')

# compare explicitly to diffusion theory by normalizing to the diffusion theory prediction
plt.figure()
for locus in range(L):
    y,x = np.histogram(allele_frequencies[:,locus], bins=af_bins, normed='True')

    # calculate the diffusion theory single locus result
    diffusion_theory = bc**(2*N*mu-1)*(1-bc)**(2*N*mu-1)*np.exp(2*N*s[locus]*bc)

    # normalize
    diffusion_theory /= np.sum(diffusion_theory)*(bc[1]-bc[0])

    # plot normalized
    plt.plot(bc, y/diffusion_theory, color=plt.cm.jet(locus*4), ls='-')

plt.title('Comparison to diffusion theory for $rN='+str(r*N)+'$, $\mu N='+str(mu*N)+'$, $N='+str(N)+'$')
plt.xlabel(r'Allele frequency $\nu$')
plt.ylabel('Allele frequency distribution/Diffusion theory result')

plt.ion()
plt.show()
