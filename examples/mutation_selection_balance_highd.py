'''
author:     Richard Neher
date:       11/07/12
content:    Example on the steady state distribution of allele frequency in a 
            balance between mutation and genetic drift using haploid_highd.
'''
# Import module (setting the path should not be necessary when the module is installed in the python path
import sys
sys.path.insert(0,'../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h

### specify parameters
N=500                               #Population size
L=64                                #number of loci
s = np.linspace(-2 ,2,L)/N          #additive selection coefficients for L loci, scaled to N
mu = 0.5/N                          #mutation rate, scaled relative to N
r  = 50.0/N/L                       #recombination rate for each interval between loci. 

### set up
pop = h.haploid_highd(L)            #produce an instance of haploid_highd with L loci
pop.carrying_capacity = N           #set the population size
pop.set_fitness_additive(0.5*s)     #set and additive fitness function. Note that FFPopSim models fitness landscape 
                                    #in a +/- rather than 0/1 basis, hence the factor 1/2
pop.recombination_model = h.CROSSOVERS                      
pop.outcrossing_rate = 1            #assign the recombination rates
pop.crossover_rate = r              #assign the recombination rates
pop.mutation_rate = mu              #assign the mutation rate

pop.set_allele_frequencies(0.5*np.ones(L),N)        #initialize the population with N individuals. 
                                                    #Population is in linkage equilibrium with frequencies 0.5
pop.evolve(10*N)                    #run for 10N generations to equilibrate



#### sample allele frequencies and store in array allele_frequencies
nsamples = 10000                    #increase for better statistics      
allele_frequencies = np.zeros((nsamples,L))

for ii in range(nsamples):
    pop.evolve(0.1*N)               #N/10 generations between successive samples
    if (ii%100==0): print ii, "out of", nsamples, "population size: ", pop.population_size, pop.number_of_clones
    allele_frequencies[ii,:]=pop.get_allele_frequencies()


### Make allele frequency histogram
af_bins = np.linspace(0,1,26)              #bins for histogram
af_bins[0]-=1.0/N
af_bins[-1]+=1.0/N
bc = 0.5*(af_bins[1:]+af_bins[:-1])    #bin centers for plotting

plt.figure()
for locus in range(L):
    y,x = np.histogram(allele_frequencies[:,locus], bins=af_bins, normed='True')
    plt.plot(bc, y, color=plt.cm.jet(locus*4))

plt.title('Comparison to diffusion theory for $rN='+str(r*N)+'$, $\mu N='+str(mu*N)+'$, $N='+str(N)+'$')
plt.text(0.3,3,"Color indicates selection coefficient \nfrom Ns=-2..2")
plt.xlabel(r'Allele frequency $\nu$')
plt.ylabel('Allele frequency distribution $f(\nu)$')


# compare explicitly to diffusion theory by normalizing to the diffusion theory prediction
plt.figure()
for locus in range(L):
    y,x = np.histogram(allele_frequencies[:,locus], bins=af_bins, normed='True')
    diffusion_theory = bc**(2*N*mu-1)*(1-bc)**(2*N*mu-1)*np.exp(2*N*s[locus]*bc)    #calculate the diffusion theory single locus result
    diffusion_theory /= np.sum(diffusion_theory)*(bc[1]-bc[0])                      #normalize
    plt.plot(bc, y/diffusion_theory, color=plt.cm.jet(locus*4), ls='-')

plt.title('Comparison to diffusion theory for $rN='+str(r*N)+'$, $\mu N='+str(mu*N)+'$, $N='+str(N)+'$')
plt.xlabel(r'Allele frequency $\nu$')
plt.ylabel('Allele frequency distribution/Diffusion theory result')
