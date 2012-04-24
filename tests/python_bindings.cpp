/*
author:     Fabio Zanini
date:       17 September 2011
content:    Extend Python with a C++ function using Richard's evolution library.
*/

/* Standard includes */
/* everything included via popgen.h */
using namespace std;


/* Custom includes */
#include "evolution.h"
#include "popgen.h"


int
check_normalization(int gts, double* frequencies)
{
 // check normalization of frequencies
 double sum = 0.0;
 for(int i=0; i < gts; i++)
  sum += frequencies[i];
 if(abs(sum - 1) > 0.01)
  return 1; 
 return 0;
}

/* HOPPING EVOLUTION */
int
evolve(double N, int generations, int gts, int loci, double* hoprates, int gts2, double* frequencies)
{
 // check array dimensions
 if ( (gts2 != gts) || (gts != (1<<loci)) )
  return 1; 
 /* check normalization */
 int flag = check_normalization(gts2, frequencies);
 if(flag)
  return flag;
 
 // set initial genotype according to *frequencies
 vector<index_value_pair> ingenotypes(gts2);
 for(int i = 0; i < gts2; i++) {
  ingenotypes[i].index = i;
  ingenotypes[i].val = frequencies[i];
 }

 // create the population
 haploid_pop_hop pop(loci,N);
 pop.init_genotypes(ingenotypes);
 pop.set_hopping_rate(hoprates);

 // evolve by hopping
 for (int i=0; i<generations; i++) {
  pop.hop();
 }
 
 // Only a signal is returned. The rest is done IN PLACE with numpy arrays and some SWIG magic!
 for (int gt=0; gt<gts; gt++)
   frequencies[gt] = pop.get_genotype_frequency(gt);

 // N.B.: the sub-arrays must NOT be released! They contain the data!
 return 0;
}



/* REAL EVOLUTION WITH POPULATION (CONSENSUS) SEQUENCING */
int
evolve_consensus(double N, int generations, double mutrate, int loci, double* fitness_additive, int gts2, double* frequencies, int reps, long seed)
{
 // check array dimensions
 if (gts2 != (1<<loci))
  return 1; 
 /* check normalization */
 int flag = check_normalization(gts2, frequencies);
 if(flag)
  return flag;

 /* Copy the initial genotypes and reset the original ones */
 double *freq0 = new double[gts2];
 for (int gt=0; gt < gts2; gt++) {
  freq0[gt] = frequencies[gt];
  frequencies[gt] = 0.0;
 }
 
 /* set initial genotype according to *frequencies */
 vector<index_value_pair> ingenotypes(gts2);
 for(int i = 0; i < gts2; i++) {
  ingenotypes[i].index = i;
  ingenotypes[i].val = freq0[i];
 }

 /* allocate random number */
 gsl_rng *rng;  // random number generator
 rng=gsl_rng_alloc(RNG);
 if(!seed)
  long seed = time (NULL) * getpid();
 gsl_rng_set (rng, seed); 

 /* MAIN LOOP */
 for(int r = 0; r < reps; r++) {

  /* generate random number */
  long rand = gsl_rng_uniform_int(rng,1000000);

  /* create the population */
  haploid_gt_dis pop(loci,N,rand);
  pop.init_genotypes(ingenotypes);
  // set mutation rates
  pop.set_mutation_rate(mutrate);

  /* set fitness */
  double fitness_additive[loci];
  for(int i = 0; i < loci; i++)
   fitness_additive[i] = 0.0;
  fitness_additive[1] = 0.2;
  pop.fitness.additive(fitness_additive,0);

  /* DEBUG: FIXME */
  /* It is somehow eaten up, but I do not know exactly how */
  printf("Input:\n");
  for(int i = 0; i < loci; i++)
   printf("%f ",fitness_additive[i]);
  printf("\n");

  printf("In population:\n");
  for(int gt = 0; gt < (1<<loci); gt++)
   printf("%f ",pop.fitness.get_coeff(gt));
  printf("\n");
 
  /* evolve without recombination */
  pop.evolve_norec(generations);
 
  /* get (big endian, integer) consensus */
  int consensus = 0;
  for(int i = 0; i< loci; i++) {
   if(pop.get_allele_frequency(i) > 0.5)
    consensus += (1<<(loci - i - 1));
  }

  /* update frequencies */
  frequencies[consensus] += 1.0;
 }

 /* normalize frequencies */
 for(int i = 0; i < gts2; i++)
  frequencies[i] /= reps;

 gsl_rng_free (rng);
 delete [] freq0;

 return 0;
}
