'''
author:     Richard Neher
date:       11/07/12
content:    compare the speed of evolution at different recombination rates
'''
# Import module
import sys
sys.path.append('../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as ffpop
import argparse

plt.ion()

#parse the command line arguments
parser = argparse.ArgumentParser(description="compare the speed of evolution at different recombination rates")
parser.add_argument('--pop', default=10000, type=float, help='Population size (N)')
parser.add_argument('--mut', default=0.0001,type=float, help='mutation rate')
parser.add_argument('--sel', default=0.01,type=float, help='effect size of beneficial mutations')
parser.add_argument('--Ttraj', default=1000,type=int, help='Length of trajectory in generations')
parser.add_argument('--dt',  default=10,type =int, help='time increments of trajectory')
params=parser.parse_args()

#set up the population
L=256                                           #the number of loci
pop=ffpop.haploid_highd(L)                      #produce an instance of haploid_highd
pop.recombination_model = ffpop.CROSSOVERS      #specify crossover recombination
pop.crossover_rate = 2.0/L                      #on average 2 crossovers per recombination event
pop.set_trait_additive(np.ones(L)*params.sel)   #set additive contribution of each locus to s
pop.mutation_rate=params.mut                    #set mutation rate
outcrossing_rates = np.linspace(0,.2, 6)        #produce an array containing the out-crossing rates to be simulated

#loop over the outcrossing rates
for r in outcrossing_rates:
    pop.outcrossing_rate=r                      #set outcrossing rates
    pop.set_wildtype(params.pop)                #initialize a population of size params.pop with all wildtype genomes
    
    #Simulate for params.Ttraj generations and record the mean fitness, fitness variance, participation ratio and the number of clones
    pfit = pop.get_fitness_statistics()
    popstat = []
    for gen in range(params.dt,params.Ttraj, params.dt):
        #append current statistics to the list
        pfit = pop.get_fitness_statistics()
        popstat.append([gen,pfit.mean, pfit.variance, pop.get_participation_ratio(), pop.number_of_clones])
        
        #evolve for dt generations and clean up
        pop.evolve(params.dt)
        pop.unique_clones()
        pop.calc_stat()
    
    #cast population statistics to an array to allow slicing
    popstat=np.array(popstat)
    
    #plot quantities of interest
    plt.figure(1)
    plt.plot(popstat[:,0],popstat[:,-2], label='r='+str(r))
    
    plt.figure(2)
    plt.plot(popstat[:,0],popstat[:,1], label='r='+str(r))
    plt.figure(3)
    plt.plot(popstat[:,0],popstat[:,2], label='r='+str(r))


#label plots and add legends
plt.figure(1)
plt.xlabel('Time')
plt.legend(loc=2)
ax=plt.gca()
ax.set_yscale('log')
plt.ylabel('Participation ratio')

plt.figure(2)
plt.legend(loc=2)
plt.ylabel('mean fitness')
plt.xlabel('Time')

plt.figure(3)
plt.legend(loc=2)
plt.ylabel('fitness variance')
plt.xlabel('Time')








