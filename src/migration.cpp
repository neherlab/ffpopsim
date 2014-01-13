#include "migration.h"
#include "ffpopsim_highd.h"

multi_population::multi_population(int loci, int two, int number_of_traits, bool three, int number_of_locations, double m_r, double m_v)
{
    //initialize all the vectors
    areal = number_of_locations;
    generation = 0;
    location.resize(number_of_locations);
    fitness_correction.resize(number_of_locations);


    //set-up the sub-populations and fitness correction
    for (int i = 0; i < number_of_locations; i++){
        fitness_correction[i] = 1.0;
        location[i] = new haploid_highd (loci, two, number_of_traits, three);
    }

    //set the migration properties
    migration_rate = m_r;
    migration_velocity = m_v;

    //set_generation
    generation = -1;

    //set initial genealogy
    genealogy.clear();

}


multi_population::multi_population(int number_of_locations, double m_r, double m_v)
{
    //initialize all the vectors
    areal = number_of_locations;
    generation = -1;
    location.resize(number_of_locations);
    fitness_correction.resize(number_of_locations);


    //set-up the sub-populations and fitness correction
    for (int i = 0; i < number_of_locations; i++){
        fitness_correction[i] = 1.0;
        location[i] = new haploid_highd (NUMBER_OF_LOCI, 0, NUMBER_OF_TRAITS, false);
    }

    //set the migration properties
    migration_rate = m_r;
    migration_velocity = m_v;


    //set_generation
    generation = -1;

    //set initial genealogy
    genealogy.clear();
}

int multi_population::population_size()
{
    int population_size = 0;
    for (int i = 0; i < areal; i ++)
    {
        population_size += point_subPopulation(i)->population_size;
    }
    return population_size;
}

multi_population::~multi_population()
{
   for (int i = 0; i < areal; i ++){
         delete location[i];
    }
    location.clear();
    fitness_correction.clear();
}


/**
 * @brief multi_population::migrate
*Migration module below
*
*
*
*
**/



void multi_population::migrate(int source)
{
    int m_n = determine_number_of_migrants(location[source]);
    if (m_n == 0) return;
    for (int i = 0; i < m_n; i ++)
    {


    //set the migrant

        int migrant_clone_No = pickup_migrant(location[source]);


    //find a detination for the migrant
    int destination;
    if (areal == 1) destination = 0;
    else destination = gsl_rng_uniform_int(location[source]->evo_generator, areal);



    //add one genotype
    location[source]->allele_frequencies_up_to_date = false;
    location[destination]->allele_frequencies_up_to_date = false;
    if (location[destination]->available_clones.size() == 0)
        location[destination]->provide_at_least(1);
    int new_gt = location[destination]->available_clones.back();
    location[destination]->available_clones.pop_back();

    location[destination]->population[new_gt].genotype = location[source]->population[migrant_clone_No].genotype;
    location[destination]->population[new_gt].clone_size = 1;

    // Update parameters for tree
    location[destination]->population[new_gt].index = new_gt;
    location[destination]->population[new_gt].location = destination;
    location[destination]->population[new_gt].parent_index = location[source]->population[migrant_clone_No].parent_index;
    location[destination]->population[new_gt].parent_location = location[source]->population[migrant_clone_No].parent_location;



    //Update fitness of the new clone
    location[destination]->calc_individual_traits(location[destination]->population[new_gt]);
    location[destination]->calc_individual_fitness_from_traits(location[destination]->population[new_gt]);
    location[destination]->check_individual_maximal_fitness(location[destination]->population[new_gt]);

    //location[destination]->population_size ++;
    location[destination]->last_clone = (new_gt < location[destination]->last_clone)?location[destination]->last_clone:new_gt;



    //Update general sub-population parameters
    location[source]->population_size --;
    location[source]->population[migrant_clone_No].clone_size -- ;
    if (location[source]->population[migrant_clone_No].clone_size < 1)
    {
        location[source]->population[migrant_clone_No].clone_size = 0;
        location[source]->available_clones.push_back(migrant_clone_No);
        location[source]->number_of_clones -- ;
    }
    location[destination]->number_of_clones ++;
    location[destination]->population_size ++;

    }
        return;
}


//choose a migrant in the population <vector> at random
int multi_population::pickup_migrant(haploid_highd* sub_population)
{
    int migrant_No = gsl_rng_uniform_int(sub_population->evo_generator, sub_population->population_size);
    int i = 0;
    int j = -1;
    while (i < migrant_No && j < sub_population->population_size)
    {
        j ++;
        i += sub_population->population[j].clone_size;

    }
    if (j < 0) j = 0;
    return j;
}





int multi_population::determine_number_of_migrants(haploid_highd* sub_population)
{
    int number_of_migrants;
    number_of_migrants = sub_population->N() * migration_rate;
    return number_of_migrants;
}


/**
 * @brief haploid_highd::evolve
 * @param gen
 * @return
 *  evolve() function is overloaded in the multi_population class
 *  Actually, this function calls the original methods from each of the sub-populations (hh)
 *
 **/



int multi_population::evolve(int loc, int gen) {
    if (HP_VERBOSE) cerr<<"multi_population::evolve(int gen)...";

    int err=0, g=0;
    location[loc]->allele_frequencies_up_to_date = false;

    /*########################################
     *
     *NO RECOMBINATION SO FAR
    // calculate an effective outcrossing rate to include the case of very rare crossover rates.
    // Since a recombination without crossovers is a waste of time, we scale down outcrossing probability
    // and scale up crossover rate so that at least one crossover is guaranteed to happen.
    if (recombination_model==CROSSOVERS)
        outcrossing_rate_effective = outcrossing_rate * (1 - exp(-number_of_loci * crossover_rate));
    else
        outcrossing_rate_effective = outcrossing_rate;

    ############################################*/

    // evolve cycle
    while((err == 0) && (g < gen)) {
        if (HP_VERBOSE) cerr<<"generation "<<generation<<endl;
        location[loc]->random_sample.clear();			//discard the old random sample
        if(err==0) err=location[loc]->select_gametes();	//select a new set of gametes (partitioned into sex and asex)
        else if(HP_VERBOSE) cerr<<"Error in select_gametes()"<<endl;
        sort(location[loc]->available_clones.begin(), location[loc]->available_clones.end(), std::greater<int>()); //sort clones in order to use the first ones again and again

        /*#########################################
         *NO RECOMBINATION SO FAR
         *if(err==0) err=add_recombinants();	//do the recombination between pairs of sex gametes
        else if(HP_VERBOSE) cerr<<"Error in recombine()"<<endl;
        ############################################3
        */

        if(err==0) err=location[loc]->mutate();		//mutation step
        else if(HP_VERBOSE) cerr<<"Error in mutate()"<<endl;
        location[loc]->random_sample.clear();			//discard the old random sample

        g++;
        //generation++;

        //add the current generation to the genealogies and prune (i.e. remove parts that do not contribute the present.


    }
    if (HP_VERBOSE) {
        if(err==0) cerr<<"done."<<endl;
        else cerr<<"error "<<err<<"."<<endl;
    }

    return err;
}




/**
 *Access functions
 *
 *
 *
 *
**/

haploid_highd* multi_population::point_subPopulation(int loc_num)
{
    return this->location[loc_num];
}


/**
 *
 *
 *Genealogy
 *
 *
 **/

//*Read all the parameters from clone and create a node_t to submit to the tree
node_t* multi_population::create_leaf_from_clone(clone_t* clone)
{
    node_t * new_leaf  = new node_t;



    //set node parameters    for the new  leaf
    new_leaf->location = clone->location;
    //new_leaf->number_of_offsprings = 0;
    new_leaf->clone_size = clone->clone_size;
    new_leaf->crossover[0] = 0;
    new_leaf->crossover[1] = 1;


    //create own_key for this node
    new_leaf->own_key.index = clone->index;
    new_leaf->own_key.location = new_leaf->location;
    new_leaf->own_key.age = generation;

    //create parent_key for this node
    new_leaf->parent_node.index = clone->parent_index;
    new_leaf->parent_node.location = clone->parent_location;
    new_leaf->parent_node.location = generation - 1;
    return new_leaf;
}



/**
 *
 *
 *
 *Genealogy
 *
 *
 ***/

///****Create new node_t for each clone in subpopulation and submit to the tree->generation_to_tree vector
void multi_population::submit_subpopulation_to_genealogy_vector(philogenetic_tree* tree, haploid_highd* sub_pop)
{
    node_t* new_leaf;

   for (int i = 0; i < sub_pop->population_size; i ++)
   {
       // check for non-zero clone sizes; all others are dummies
       if (sub_pop->population[i].clone_size != 0)
       {
           new_leaf = create_leaf_from_clone(&sub_pop->population[i]);
           tree->push_node_to_generation_vector(new_leaf);

       }
   }
   return;
}


//***Submit each sub-population to the tree->generation_to_tree vector
void multi_population::submit_population_to_genealogy_vector(philogenetic_tree* tree)
{
    for (int i = 0; i < areal; i ++)
    {
        submit_subpopulation_to_genealogy_vector(tree, point_subPopulation(i));
    }

}
/*
//***Submit each sub-population to the tree->generation_to_tree vector
void multi_population::submit_vector_to_tree (int locus)
{
    if (locus < 0 or locus > genealogy.size())
        return;
    genealogy[locus].    add_generation(, baseline)
}

*/

int multi_population::track_locus_genealogy(vector < int > gen_loci)
{
    genealogy.resize(gen_loci.size());

    //Note: you must track genealogies BEFORE the population is set
    if((generation != -1)){

        return -1;
    }
    for (int i = 0; i < areal; i ++)
    {
        if (point_subPopulation(i)->get_number_of_clones() != 0)
            return -1;
    }

    track_genealogy=true;

    //genealogy.reset();
    /*for (unsigned int i=0; i<gen_loci.size(); i++){
        genealogy[i].track_locus(gen_loci[i]);
        genealogy[i].extend_storage(population_size());
    }*/

    if (HP_VERBOSE){cerr<<"done\n";}
    return 0;
}


int multi_population::get_tree_size(int tree_num)
{
    if (tree_num > 0 and tree_num < genealogy.size())
    {

        return 0;
    }
}


string multi_population::print_newick_genealogy(int genlocus)
{
    if (genlocus > 0 and genlocus < genealogy.size())
    {
        return genealogy[genlocus].print_newick();
    }
    return "\n";
}


philogenetic_tree * multi_population::point_tree(int genlocus)
{
    if ((genlocus > -1) and (genlocus < genealogy.size()))
    {
        return &genealogy[genlocus];
    }
    return NULL;
}

