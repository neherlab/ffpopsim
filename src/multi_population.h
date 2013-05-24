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
    //vector <haploid_highd> sub_population;
    haploid_highd * sub_population;
    //vector <haploid_highd> sub_population;
    int population_size();

    int number_of_locations;
    int track_genealogy;
    double fitness_max;
    int generation;

    double migration_rate;

    void add_migrating_clone_to_genealogy(int locusIndex, int old_location, int new_location,  int dest, int parent, int left, int right, int cs, int n);
    int mutate(int location);
    unsigned int flip_single_locus(int location, unsigned int clonenum, int locus);
    unsigned int flip_single_locus(int location, int locus);




public:

    int get_locations () {return number_of_locations;};

    multi_population(int new_locations, int L_in, int n_o_traits = 1);
    ~multi_population();

    int set_theonly_wildtype(int new_location, int new_N);

    multi_locus_genealogy genealogy;
    int track_locus_genealogy(vector<int> loci);
    int submit_subpop_genealogy(int sub_pop_No);
    int submit_pop_genealogy();
    double max_fitness();

    int migrate();
    int migrate(int source);
    int determine_number_of_migrants(haploid_highd sub_population);
    int pickup_migrant(haploid_highd sub_population);
    int transfer_clone(int sub_pop_source, int sub_pop_destination, int source);

    int number_of_migration_events;

    int determine_number_of_migrants(int sub_pop_No);
    int determine_migration_destination();
    int determine_migrant(int sub_pop_num);



    void reset();
    haploid_highd * point_sub_pop(int i){return &sub_population[i];};
    void set_global_generation(int generation);


    int evolve(int location, int gen);

    int set_migration_rate(double new_rate){migration_rate = new_rate; return 0;};


};


#endif // MULTI_POPULATION_H
