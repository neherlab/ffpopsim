# Import module
import sys
sys.path.append('../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as ffpop
import argparse
import pickle

#parse the command line arguments
parser = argparse.ArgumentParser(description="Simulate a population on a mixed additive/epistatic fitness function")
parser.add_argument('--pop', default=10000, type=float, help='Population size (N)')
parser.add_argument('--rec',  default=0,type=float, help='out-crossing rate (r)')
parser.add_argument('--sigma', default=0.05,type=float, help='Sigma')
parser.add_argument('--hsq', default=0,type=float, help='heritability')
parser.add_argument('--Ttraj', default=200,type=float, help='Length of trajectory in generations')
parser.add_argument('--dt',  default=1,type =int, help='time increments of trajectory')
parser.add_argument('--outdir',  default='./', help='Directory into which output is directed')
parser.add_argument('--runs',  default=1, type = int , help='number of runs')
params=parser.parse_args()

#set up the population
L=64
pop=ffpop.haploid_highd(L)
pop.outcrossing_rate=params.rec
pop.set_random_epistasis(params.sigma*np.sqrt(1-params.hsq))
pop.recombination_model = ffpop.FREE_RECOMBINATION
if (params.hsq>0):
    pop.set_trait_additive(np.ones(L)*params.sigma*sqrt(params.hsp/L))


#evolve
popstat = []
#loop over the requested number of runs
for ri in xrange(params.runs):
    #initialize
    pop.set_allele_frequencies(np.ones(L)*0.5, params.pop)
    pfit = pop.get_fitness_statistics()
    popstat.append([[0,pfit.mean, pfit.variance, pop.participation_ratio, pop.number_of_clones]])
    for gen in range(0,params.Ttraj, params.dt):
        pop.evolve(params.dt)
        pop.unique_clones()
        pfit = pop.get_fitness_statistics()
        popstat[-1].append([0,pfit.mean, pfit.variance, pop.participation_ratio, pop.number_of_clones])
    

popstat=np.array(popstat)
pickle_file = open(params.outdir+'/statistics_r_'+str(params.rec)+'_N_'+str(params.pop)+'.pickle', 'w')
pickle.dump((params, popstat), pickle_file)
pickle_file.close()




