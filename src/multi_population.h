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

/**
 * @brief Population class for high-dimensional simulations of multiple migrating populations in changing environmental condditions.
 *
 * This class is the main object for simulating multiple geographically divided population in changing environment.
 * Each population is modelled by a single haploid_highd class enabling manipulation of the populations with long genomes (\f$L\f$) larger than 20.
 *
 */
class multi_population
{
<<<<<<< Updated upstream
    double MAX_MIGRATION_RATE;
    //vector <haploid_highd> sub_population;
    vector <haploid_highd*> sub_population;
    //vector <haploid_highd> sub_population;
    int L;
=======
    vector <haploid_highd*> sub_population; //main storage for sub_populations
>>>>>>> Stashed changes
    int number_of_locations;
    int number_of_loci;
    int track_genealogy; //0=no genealogy; 1=only local trees are stored; 2=global trees.
    double fitness_max;//maximal fitness among all sub-populations
    double migration_rate;//global migration rate

    void add_migrating_clone_to_genealogy(int locusIndex, int old_location, int new_location,  int dest, int parent, int left, int right, int cs, int n);
    int mutate(int location);
    unsigned int flip_single_locus(int location, unsigned int clonenum, int locus);
    unsigned int flip_single_locus(int location, int locus);
    int generation;
    int determine_number_of_migrants(int sub_pop_No);
    int determine_migration_destination();
    int determine_migrant(int sub_pop_num);

<<<<<<< Updated upstream
    //genealogy
    multi_locus_genealogy genealogy;
    int submit_subpop_genealogy(int sub_pop_No);
    int submit_pop_genealogy();
    void set_global_generation(int generation);

    //evolution
    int evolve_local(int location, int gen=1);




public:

    void reset();
    multi_population(int new_locations, int L_in, int n_o_traits = 1, int rng_seed = 0);
    ~multi_population();
    haploid_highd * point_sub_pop(int i){return sub_population[i];};

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
    int track_locus_genealogy(vector<int> loci);
    string print_newick();


    //evolution
    int evolve(int gen = 1);
    int migrate();
    int migrate(int source);


    //init population
    void add_random_genotype(int N_in);
=======
    int determine_number_of_migrants(int sub_pop_No); //number of migrants in the sub_pop_No location. Determined using the migration rate.
    int determine_migration_destination();//Find random destination for a migration event
    int determine_migrant(int sub_pop_num);//pick a random clone to migrate
    int submit_subpop_genealogy(int sub_pop_No);//submits sub-population to the local genealogy
    int evolve_local(int location, int gen); //main function for sub-population evolution.
    int migrate();//migrate from randomly determined source to randomly determined destination


public:
//temporarily public stuff
    int generation;
    void set_global_generation(int generation);//controls generation counters in sub-populations
    int submit_pop_genealogy(); //submit all local genealogies to the global genealogy
    int migrate(int source); //migrate from specified source to randomly determined destination


    multi_population(int new_locations, int L_in, int n_o_traits = 1, int rng_seed = 0);
    ~multi_population();
    int get_locations () {return number_of_locations;};
    int track_locus_genealogy(vector<int> loci); //enables genealogy
    double max_fitness();//get maximal fitness
    multi_locus_genealogy genealogy;//class for storage the global genealogy data. Collects all the data from local genealogies stored in haploid_highd


    int evolve(int gen = 1);
    int number_of_migration_events; //migration events counter
    void reset();
    haploid_highd * point_sub_pop(int i){return sub_population[i];};//returns pointer to a haploid-highd sub-population

    int get_global_generation(){return generation;};
    int set_migration_rate(double new_rate){migration_rate = new_rate; return 0;};
    double get_migration_rate(){return migration_rate;};

    void set_carrying_capacity(int N);
    int get_carrying_capacity(){return point_sub_pop(0)->carrying_capacity;};
    int population_size();

    void set_mutation_rate(double mu);
    double get_mutation_rate(){return point_sub_pop(0)->get_mutation_rate();};

    void set_outcrossing_rate(double o_rate);
    double get_outcrossing_rate(){return point_sub_pop(0)->outcrossing_rate;};

    void set_crossover_rate(double c_rate);
    double get_crossover_rate(){return point_sub_pop(0)->crossover_rate;};

    void set_recombination_model(int rec_model);
    int get_recombination_model(){return point_sub_pop(0)->recombination_model;};

    void set_trait_coefficients(int coefficient, vector <int> loci, int trait_no);
    void set_trait_weights(double* new_weights);

    void update_traits();
    void update_fitness();
    int add_random_genotype(int N=1);


>>>>>>> Stashed changes

protected:

    int determine_number_of_migrants(haploid_highd sub_population);
    int pickup_migrant(haploid_highd sub_population);
    int transfer_clone(int sub_pop_source, int sub_pop_destination, int source);




};


#endif // MULTI_POPULATION_H
