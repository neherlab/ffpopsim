#ifndef MULTI_POPULATION_H
#define MULTI_POPULATION_H
#include "ffpopsim_highd.h"
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <boost/dynamic_bitset.hpp>

class multi_population
{
    double MAX_MIGRATION_RATE;
    //vector <haploid_highd> sub_population;
    vector <haploid_highd*> sub_population;
    //vector <haploid_highd> sub_population;
    int L;
    int number_of_locations;
    int track_genealogy;
    double fitness_max;

    double migration_rate;

    void add_migrating_clone_to_genealogy(int locusIndex, int old_location, int new_location,  int dest, int parent, int left, int right, int cs, int n);
    int mutate(int location);
    unsigned int flip_single_locus(int location, unsigned int clonenum, int locus);
    unsigned int flip_single_locus(int location, int locus);
    int generation;
    int determine_number_of_migrants(int sub_pop_No);
    int determine_migration_destination();
    int determine_migrant(int sub_pop_num);

    int submit_subpop_genealogy(int sub_pop_No);
    int submit_pop_genealogy();
    void set_global_generation(int generation);
    int evolve_local(int location, int gen);
    void migrate();
    void migrate(int source);

public:


    void reset();
    haploid_highd * point_sub_pop(int i){return sub_population[i];};
    multi_population(int new_locations, int L_in, int n_o_traits = 1, int rng_seed = 0);
    ~multi_population();

    //parameteres
    int get_locations () {return number_of_locations;};
    int N();
    int get_generation(){return generation;}
    double max_fitness();
    int number_of_migration_events;

    void set_migration_rate(double new_rate){migration_rate = new_rate;};
    double get_migration_rate(){return migration_rate;};

    void set_mutation_rate(double mu);
    double get_mutation_rate(){return point_sub_pop(0)->get_mutation_rate();};

    void set_carrying_capacity(int capacity);
    int get_carrying_capacity(){return point_sub_pop(0)->carrying_capacity;};

    void set_outcrossing_rate(double o_rate);
    double get_outcrossing_rate(){return point_sub_pop(0)->outcrossing_rate;}

    void set_crossover_rate(double c_rate);
    double get_crossover_rate(){return point_sub_pop(0)->crossover_rate;};

    void set_recombination_model(int r_model);
    int  get_recombination_model(){return point_sub_pop(0)->recombination_model;};

    void set_trait_coefficient(double coefficient, vector<int> loci, int trait_no);
    void set_trait_weights(double* weights);

    //genealogy
    multi_locus_genealogy genealogy;
    void track_locus_genealogy(vector<int> loci);




    //evolution

    int evolve(int gen = 1);


    //init population
    void set_random_genotype(int N_in);

protected:

    int determine_number_of_migrants(haploid_highd sub_population);
    int pickup_migrant(haploid_highd sub_population);
    int transfer_clone(int sub_pop_source, int sub_pop_destination, int source);




};


#endif // MULTI_POPULATION_H
