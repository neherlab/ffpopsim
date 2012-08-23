import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h

### specify parameters
N=500                               #Population size
L=4                                 #number of loci
s = np.linspace(-0.2 ,0.7,L)/N      #additive selection coefficients for L loci, scaled to N
mu = 0.4/N                          #mutation rate, scaled relative to N
r = 5.0/N*np.ones(L-1)              #recombination rate for each interval between loci. 

print "N =",N, ', L =',L, ', mu N = ', mu*N, ', s N = ', N*s 

### set up
pop = h.haploid_lowd(L)             #produce an instance of haploid_lowd with L loci
pop.carrying_capacity = N           #set the population size
pop.set_fitness_additive(0.5*s)     #set and additive fitness function. Note that FFPopSim models fitness landscape 
                                    #in a +/- rather than 0/1 basis, hence the factor 1/2
pop.set_recombination_rates(r)      #assign the recombination rates
pop.set_mutation_rates(mu)           #assign the mutation rate

pop.set_genotypes([0],[N])          #initialize the population with N individuals with genotypes 0, that is ----
pop.evolve(10*N)                    #run for 10N generations to equilibrate



#### sample allele frequencies and store in array allele_frequencies
nsamples = 10000                    #increase for better statistics      
allele_frequencies = np.zeros((nsamples,L))

for ii in range(nsamples):
    pop.evolve(0.1*N)               #N/10 generations between successive samples
    allele_frequencies[ii,:]=pop.get_allele_frequencies()


### Make allele frequency histogram
af_bins = np.linspace(0,1,26)                   #bins for histogram
bin_centers = 0.5*(af_bins[1:]+af_bins[:-1])    #bin centers for plotting
nu = np.linspace(0.01,0.99,99)                  #denser frequency grid for theory curves

cols = ['r', 'b', 'g', 'm', 'c']
for locus in range(L):
    y,x = np.histogram(allele_frequencies[:,locus], bins=af_bins, normed='True')
    plt.plot(bin_centers, y, color=cols[locus], label = r'$N s_'+str(locus+1)+'='+str(N*s[locus])+'$')
    
    diffusion_theory = nu**(2*N*mu-1)*(1-nu)**(2*N*mu-1)*np.exp(2*N*s[locus]*nu)    #calculate the diffusion theory single locus result
    diffusion_theory /= np.sum(diffusion_theory)*(nu[1]-nu[0])                      #normalize
    plt.plot(nu,diffusion_theory, color=cols[locus], ls='--')                       #add to plot with same color and dashed lines

plt.legend(loc=9)
plt.title('Comparison to diffusion theory for $rN='+str(r[0]*N)+'$, $\mu N='+str(mu*N)+'$, $N='+str(N)+'$')
plt.xlabel(r'Allele frequency $\nu$')
plt.ylabel('Allele frequency distribution')
