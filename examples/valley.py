# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       25/05/12
content:    Find the time for valley crossing
'''
# Import module
import sys
sys.path.append('../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h


# Globals
simulate = False                # Use extant data or repeat the simulation?
L = 4                           # Number of loci
N = 1e10                        # Population size
s1 = 1e-5                       # Fitness advantage of wildtype
s2 = 0.01                       # Fitness advantage of quadruple mutant

# Recombination rates to check out
rs = np.logspace(-4,-3,10).tolist() + \
     [0.00125, 0.0015, 0.00175, 0.002,\
      0.00225, 0.0025, 0.00275, 0.003,\
      0.0031, 0.0032, 0.0033, 0.0034, \
      0.0035, 0.00375, 0.004, 0.005, \
      0.0075, 0.01]

# Mutation rates to check out
mutation_rates=[1e-7,1e-6, 1e-5]



# Script
if __name__ == '__main__':

    # Prepare figure
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    lstyles = [[1,0],
               [2, 4, 5, 4],
               [7, 4, 7, 4]]
    colors=['b', 'g', 'r', 'm', 'c']
    
    # Cycle over mutation rates
    for k, mu in enumerate(mutation_rates):

        # Simulate the population? ...
        if simulate:

            # Prepare result vectors
            times = np.zeros_like(rs)
            dtimes = np.zeros_like(rs)

            for i, r in enumerate(rs):    
                ttmp = []
                for j in xrange(50):
        
                    c = h.haploid_lowd(L)
                    c.set_genotypes([0],[N])
                    c.set_recombination_rates(r*np.ones(c.L-1))
                    c.set_mutation_rates(mu)
                    c.set_fitness_function([0b0, 0b1111], [s1, s1+s2])
                    c.evolve()
                    # cross valley
                    gens = 100
                    while c.get_genotype_frequency(0b1111)<0.5 and c.generation<1e6:
                        c.evolve(gens)
                    
                    if (c.generation<1e6):
            	       	ttmp.append(c.generation)
                    else: 
                        ttmp.append(np.nan)
                        break
                    
                times[i] = np.mean(ttmp)
                dtimes[i] = np.std(ttmp)
                print 'r = '+str(r)+'\tTime to cross the valley: '+str(times[i])+' generations' 
        
            # Save results in a text file
            np.savetxt('valley_mu_'+str(mu)+'.dat',zip(rs,times, dtimes))

        # or read results from file?
        else:
            rs, times, dtimes = np.loadtxt('valley_mu_'+str(mu)+'.dat', unpack=True)

        # Plot
        line = ax.errorbar(rs, times, dtimes,
                    c=colors[k],
                    lw=2,
                    label=r'$\mu=10^{'+str(int(np.log10(mu)))+'}$')
        line[0].set_dashes(lstyles[k])

    # Finish image
    ax.set_title(r'Population size $N=10^{'+str(int(np.log10(N)))+'}$')
    ax.set_xlabel(r'r [gen$^{-1}$]')
    ax.set_ylabel(r'crossing time [gen]')
    ax.set_title('Time for valley crossing')
    ax.legend(loc=9)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Show the plot
    plt.ion()
    fig.show()
