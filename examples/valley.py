# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       25/05/12
content:    Find the time for valley crossing
'''
# Import modules (setting the path should not be necessary when the module is
# installed in the PYTHONPATH)
import sys
sys.path.insert(0, '../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h


# specify parameters
simulate = False                # use extant data or repeat the simulation?
L = 4                           # number of loci
N = 1e10                        # population size
s1 = 1e-5                       # fitness advantage of wildtype (half)
s2 = 0.01                       # fitness advantage of quadruple mutant (half)

# recombination rates to check out
rs = np.logspace(-4,-3,10).tolist() + \
     [0.00125, 0.0015, 0.00175, 0.002,\
      0.00225, 0.0025, 0.00275, 0.003,\
      0.0031, 0.0032, 0.0033, 0.0034, \
      0.0035, 0.00375, 0.004, 0.005, \
      0.0075, 0.01]

# mutation rates to check out
mutation_rates=[1e-7, 1e-6, 1e-5]


# script
if __name__ == '__main__':

    # prepare figure
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    lstyles = [[1,0],
               [2, 4, 5, 4],
               [7, 4, 7, 4]]
    colors=['b', 'g', 'r', 'm', 'c']
    
    # cycle over mutation rates
    for k, mu in enumerate(mutation_rates):

        # simulate the population? ...
        if simulate:

            # prepare result vectors
            times = np.zeros_like(rs)
            dtimes = np.zeros_like(rs)

            for i, r in enumerate(rs):    
                ttmp = []
                for j in xrange(50):
        
                    # set up population
                    pop = h.haploid_lowd(L)         # produce an instance of haploid_lowd with L loci

                    pop.set_genotypes([0], [N])     # set a wildtype population of size N

                    pop.set_recombination_rates(r)  # recombination rate
                    pop.set_mutation_rates(mu)      # mutation rate

                    # set the fitness valley:
                    #    - intermediates have fitness 0
                    #    - wildtype has fitness 2 s1
                    #    - quadruple mutant has fitness 2 (s1 + s2)
                    # Note: FFPopSim reasons in the +/- basis.
                    pop.set_fitness_function([0b0, 0b1111],
                                             [s1, s1+s2])

                    # cross valley
                    gens = 100
                    while (pop.get_genotype_frequency(0b1111) < 0.5) and (pop.generation < 1e6):
                        pop.evolve(gens)
                    
                    if (pop.generation < 1e6):
            	       	ttmp.append(pop.generation)
                    else: 
                        ttmp.append(np.nan)
                        break
                    
                # store crossing times
                times[i] = np.mean(ttmp)
                dtimes[i] = np.std(ttmp)
                print 'r = '+str(r)+'\tTime to cross the valley: '+str(times[i])+' generations' 
        
            # save results in a text file
            np.savetxt('valley_mu_'+str(mu)+'.dat',zip(rs,times, dtimes))

        # ...or read results from file?
        else:
            rs, times, dtimes = np.loadtxt('valley_mu_'+str(mu)+'.dat', unpack=True)

        # plot
        line = ax.errorbar(rs, times, dtimes,
                           c=colors[k],
                           lw=2,
                           label=r'$\mu=10^{'+str(int(np.log10(mu)))+'}$')
        line[0].set_dashes(lstyles[k])

    ax.set_title(r'Population size $N=10^{'+str(int(np.log10(N)))+'}$')
    ax.set_xlabel(r'r [gen$^{-1}$]')
    ax.set_ylabel(r'crossing time [gen]')
    ax.set_title('Time for valley crossing')
    ax.legend(loc=9)
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.ion()
    fig.show()
