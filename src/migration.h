#ifndef MIGRATION_H
#define MIGRATION_H
#include "ffpopsim_highd.h"
#include "philogenetic_tree.h"


#define NUMBER_OF_LOCI 1000
#define NUMBER_OF_TRAITS 1


class multi_population
{
    //vector of available locations
    vector<haploid_highd *> location;
    int population_size();

    //number of locations
    int areal;

    //migration function and properties
    double migration_rate;
    double migration_velocity;

    //fitness correction
    vector<double> fitness_correction;

    //Genealogy
    vector < philogenetic_tree > genealogy; //Size of the vector is the number of tracked loci

    node_t * create_leaf_from_clone(clone_t* clone);
    void submit_subpopulation_to_genealogy_vector(philogenetic_tree* tree, haploid_highd* sub_pop);
    bool track_genealogy;


public:
    //current generation
    int generation;
    int evolve(int location_num, int gen=1);


    /**
    *   Temporary public functions
    **/
    void  migrate(int source);
    int determine_number_of_migrants(haploid_highd* sub_population);
    int pickup_migrant(haploid_highd* sub_population);


    //constructor and destructor
    multi_population(int number_of_locations, double migration_rate = 1, double migration_velocity = 1.0);
    multi_population(int loci = 1000, int two = 0, int number_of_traits = 1, bool three = false, int number_of_locations =1 , double migration_rate = 1, double migration_velocity = 1.0);

    ~multi_population();


    /**
     *
     *  Access functions
     *
     **/
    haploid_highd* point_subPopulation(int location_num);



    /**
     *Genealogy
     *
     **/
    int track_locus_genealogy(vector < int > gen_loci);
    void submit_population_to_genealogy_vector(philogenetic_tree* tree);
    int get_tree_size(int tree_num);
    string print_newick_genealogy(int genlocus);

    philogenetic_tree * point_tree(int genlocus);


};






#endif // MIGRATION_H
