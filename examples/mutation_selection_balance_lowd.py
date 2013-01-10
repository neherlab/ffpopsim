'''
author:     Richard Neher, Fabio Zanini
date:       11/07/12
content:    Example on the steady state distribution of allele frequency in a 
            balance between mutation and genetic drift using haploid_lowd.
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
L = 4                               # number of loci
s = np.linspace(-0.2 ,0.7, L) / N   # additive selection coefficients for L loci, scaled to N
mu = 0.4 / N                        # mutation rate, scaled to N
r = 5.0 / N                         # recombination rate for each interval between loci. 

# set up population
pop = h.haploid_lowd(L)             # produce an instance of haploid_lowd with L loci

# set and additive fitness function. Note that FFPopSim models fitness landscape
# in a +/- rather than 0/1 basis, hence the factor 1/2
pop.set_fitness_additive(0.5 * s) 

pop.set_mutation_rates(mu)          # mutation rate
pop.set_recombination_rates(r)      # recombination rate (CROSSOVERS model by default)

# initialize the population with N wildtype individuals, that is ----
pop.carrying_capacity = N           # set the population size
pop.set_genotypes([0], [N])

print "Evolve for >> N generations and compare allele frequency distributions \nto expectations from diffusion theory."
pop.status()
pop.evolve(10 * N)                  # run for 10N generations to equilibrate

# evolve and sample allele frequencies
nsamples = 10000      
allele_frequencies = np.zeros((nsamples,L))
for ii in range(nsamples):
    pop.evolve(0.1 * N)             # N / 10 generations between successive samples

    # get allele frequencies
    allele_frequencies[ii,:] = pop.get_allele_frequencies()


# prepare allele frequency histogram
af_bins = np.linspace(0,1,26)                   # bins for histogram
bin_centers = 0.5*(af_bins[1:]+af_bins[:-1])    # bin centers for plotting
nu = np.linspace(0.01,0.99,99)                  # denser frequency grid for theory curves

# plot results
cols = ['r', 'b', 'g', 'm', 'c']
for locus in range(L):

    # make histogram
    y,x = np.histogram(allele_frequencies[:,locus], bins=af_bins, density='True')

    # plot
    plt.plot(bin_centers, y, color=cols[locus], label = r'$N s_'+str(locus+1)+'='+str(N*s[locus])+'$')
    
    # calculate the diffusion theory single locus result
    diffusion_theory = nu**(2*N*mu-1)*(1-nu)**(2*N*mu-1)*np.exp(2*N*s[locus]*nu)
    # normalize
    diffusion_theory /= np.sum(diffusion_theory)*(nu[1]-nu[0])

    # plot diffusion theory with dashed lines
    plt.plot(nu,diffusion_theory, color=cols[locus], ls='--')

plt.legend(loc=9)
plt.title('Comparison to diffusion theory for $rN='+str(r*N)+'$, $\mu N='+str(mu*N)+'$, $N='+str(N)+'$')
plt.xlabel(r'Allele frequency $\nu$')
plt.ylabel('Allele frequency distribution')

plt.ion()
plt.show()
