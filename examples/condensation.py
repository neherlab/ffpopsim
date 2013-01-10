'''
author:     Richard Neher, Fabio Zanini
date:       11/07/12
content:    Genotype condensation is driven by epistasis.
'''
# Import modules (setting the path should not be necessary when the module is
# installed in the PYTHONPATH)
import sys
sys.path.insert(0, '../pkg/python')

import time
import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as ffpop


# parse the command line arguments
import argparse
parser = argparse.ArgumentParser(description="Simulate a population on a mixed additive/epistatic fitness function")
parser.add_argument('--pop', default=10000, type=float, help='Population size (N)')
parser.add_argument('--rec',  default=0,type=float, help='out-crossing rate (r)')
parser.add_argument('--sigma', default=0.05,type=float, help='Sigma')
parser.add_argument('--hsq', default=0,type=float, help='heritability')
parser.add_argument('--Ttraj', default=200,type=int, help='Length of trajectory in generations')
parser.add_argument('--dt',  default=1,type =int, help='time increments of trajectory')
params=parser.parse_args()

# set up the population
L = 1000                                            # number of loci
pop = ffpop.haploid_highd(L)                        # create an instance of the class
pop.outcrossing_rate = params.rec                   # outcrossing rate
pop.recombination_model = ffpop.CROSSOVERS          # recombination model

# set random epistatic fitness lanscape
pop.set_random_epistasis(params.sigma * np.sqrt(1 - params.hsq))

# set heritability
if (params.hsq > 0):
    pop.set_trait_additive(np.ones(L) * params.sigma * sqrt(params.hsp / L))


# initialize population in linkage equilibrium with all
# allele frequencies at 0.5
pop.set_allele_frequencies(np.ones(L) * 0.5, params.pop)

# evolve and store statistics along the way
popstat = []
t1 = time.time()
print "\nEvolve a population to illustrate the effect of epistasis on population structure"
pop.status()
print
for gen in range(params.dt,params.Ttraj, params.dt):
    if gen%(params.Ttraj/10)==0: print gen, "out of", params.Ttraj, "generations"
    # append current statistics to the list
    pfit = pop.get_fitness_statistics()
    popstat.append([gen,pfit.mean,
                    pfit.variance,
                    pop.participation_ratio,
                    pop.number_of_clones])
    
    # evolve for dt generations and clean up
    pop.evolve(params.dt)
    pop.unique_clones()
    pop.calc_stat()

popstat = np.array(popstat)
print "Elapsed time:", round(time.time() - t1,2), 's'

# plot quantities of interest
plt.figure(1)
plt.plot(popstat[:,0], popstat[:,-2])
plt.xlabel('Time')
plt.ylabel('Participation ratio')

plt.figure(2)
plt.plot(popstat[:,0], popstat[:,-1])
plt.yscale('log')
plt.xlabel('Time')
plt.ylabel('Number of clones')

plt.figure(3)
plt.plot(popstat[:,0], popstat[:,1], label='fitness mean')
plt.plot(popstat[:,0], np.sqrt(popstat[:,2]), label='fitness standard deviation')
plt.legend(loc=2)
plt.xlabel('Time')

plt.ion()
plt.show()
