#include "multi_population.h"


class hh_source : public haploid_highd{
    double sigma;


    virtual double trait_function(double aTrait1, double aTrait2)
    {
        double t1 = aTrait1;
        double t2 = aTrait2;
        double t1_0 = cos(phi_0) * rho_0;
        double t2_0 = sin(phi_0) * rho_0;

        return  trait_weights[0] / 5 * exp(-  ( (pow(t1, 2)) + pow(t2, 2) )  / (200));
    }

public:

    //get/set methods for fitness map parameters
    void set_phi_0 (double aPhi_0){phi_0 = aPhi_0; return;}
    double get_phi_0 (){return phi_0;}
    void set_offset(double aRho){rho_0 = aRho; return;}
    double get_offset(){return rho_0;}

    //This method doeas the job of trait -> fitness conversion
    virtual void calc_individual_fitness_from_traits(clone_t &tempgt) {tempgt.fitness = trait_function(tempgt.trait[0], tempgt.trait[1]);}

    hh_source()
    {
        phi_0 = 0;
        rho_0 = 0;
        sigma = 100;
    }
    ~hh_source(){}
};

class hh_1 : public haploid_highd{
    double sigma;
    virtual double trait_function(double aTrait1, double aTrait2)
    {
        double t1 = aTrait1;
        double t2 = aTrait2;
        double rho = 0;
        //double t1_0 = cos(phi_0) * rho_0;
        //double t2_0 = sin(phi_0) * rho_0;
        rho = sqrt(pow((t1 + 1e-5), 2) + pow((t2 + 1e-5), 2));
        double d_phi = acos( ((t1 + 1e-5) * cos(phi_0) + (t2 + 1e-5) * sin(phi_0)) / rho);
        return   trait_weights[0] / 5 * exp(1) *  rho / rho_0 * exp(- rho / rho_0) * exp(- pow((d_phi), 2) / sigma);
    }
public:
    //get/set methods for fitness map parameters
    void set_phi_0 (double aPhi_0)
    {
        double new_phi = aPhi_0;
        while (new_phi < 0 )
        {
                new_phi = new_phi + 2 * 3.14;
        }

        while (new_phi > 2 * 3.14)
        {
                new_phi = new_phi - 2 * 3.14;
        }

        phi_0 = new_phi;
        return;
    }
    double get_phi_0 (){return phi_0;}
    void set_offset(double aRho){rho_0 = aRho; return;}
    double get_offset(){return rho_0;}

    //This method doeas the job of trait -> fitness conversion
    virtual void calc_individual_fitness_from_traits(clone_t &tempgt) {tempgt.fitness = trait_function(tempgt.trait[0], tempgt.trait[1]);}

    hh_1()
    {
        phi_0 = 3.14 / 4;
        rho_0 = 8;
        sigma = 0.1;
    }


    ~hh_1(){}
};



/**
  * @brief Default constructor for simulating multiple populations
  * @param new_locations number of locations to be simulated
  * @param L_in number of loci
  * @param n_o_traits number of traits
  */
multi_population::multi_population(int new_locations, int L_in, int n_o_traits, int rng_seed)
{
<<<<<<< Updated upstream
    if (new_locations < 0)
    {
        if (HP_VERBOSE)
            cerr << "too few locations provided for the required type of simulations!"<< endl;
        throw(HP_BADARG);
    }

    try {
        haploid_highd* loc;
        for (int i = 0; i < new_locations; i ++)

        {
           loc = new haploid_highd();
           loc->set_up(L_in, rng_seed, n_o_traits);
           sub_population.push_back(loc);
        }

       L = L_in;
       generation = 0;
       number_of_locations = new_locations;
       track_genealogy = 0; //No genealogy by default
       migration_rate = 0.0; //No migration by default
       number_of_migration_events = 0;
    }catch (int err)
    {
        throw err;
    }
=======


    //Location number check:
    if (new_locations < 1)
    {
        if (HP_VERBOSE)
            cerr << "You must specify at least one location!"<< endl;
        throw(HP_BADARG);
    }

   haploid_highd* loc;
   //FIXME: set the fitness landscape from outside (may be, provide uilities)
   try{
       loc = new hh_source();
       loc->set_up(L_in, rng_seed, n_o_traits);
       sub_population.push_back(loc);

       for (int i = 1; i < new_locations; i ++)
       {
           loc = new hh_1();
           loc->set_up(L_in, rng_seed, n_o_traits);
           sub_population.push_back(loc);
       }
   }catch (int e){
       throw e;
   }

    set_global_generation(0);
    number_of_locations = new_locations;
    track_genealogy = 0; //No genealogy by default
    migration_rate = 0.0; //No migration by default
    number_of_migration_events = 0;
    number_of_loci = L_in;
>>>>>>> Stashed changes
}



/**
  * @brief Memory release here
*/
multi_population::~multi_population()
{
    reset();
    for (int i = 0; i < number_of_locations; i ++)
     {
        //sub_population[i].~haploid_highd();
        //delete sub_population[i];
     }
}



<<<<<<< Updated upstream

int multi_population::N()
=======
/**
  * @brief Returns overall size of the population
  @return size of the population
  */
int multi_population::population_size()
>>>>>>> Stashed changes
{
    int size = 0;
    for (int i = 0; i < number_of_locations; i ++)
    {
        size += sub_population[i]->N();
    }

    return size;
}


void multi_population::reset()
{
    number_of_locations = 0;
    sub_population.clear();
 }

void multi_population::set_carrying_capacity(int N) {
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->carrying_capacity = N;
    }
}
void multi_population::set_mutation_rate(double mu) {
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->set_mutation_rate(mu);
    }
}
void multi_population::set_outcrossing_rate(double o_rate) {
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->outcrossing_rate = o_rate;
    }
}
void multi_population::set_crossover_rate(double c_rate) {
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->crossover_rate = c_rate;//1e-2;
    }
}

void multi_population::set_recombination_model(int rec_model){
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->recombination_model = rec_model;
    }
}

void multi_population::set_trait_coefficients(int coefficient, vector <int> loci, int trait_no)
{
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->add_trait_coefficient(0, loci, 1);
    }
}

void multi_population::set_trait_weights(double* new_weights)
{
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->set_trait_weights(new_weights);
    }
}

void multi_population::update_traits()
{
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->update_traits();
    }
}
void multi_population::update_fitness()
{
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->update_fitness();
    }
}


int multi_population::add_random_genotype(int N)
{
    if(N == 0) {
        if (HP_VERBOSE) cerr<<"the desired population size must be at least 1."<<endl;
        return HP_BADARG;
    }
    // set the carrying capacity if unset
    if(get_carrying_capacity() < HP_NOTHING)
            set_carrying_capacity(N * 100);
    //point_sub_pop(0)->allele_frequencies_up_to_date = false;

    boost::dynamic_bitset<> genotype(number_of_loci);
    bool r = 0;
    for( int i = 0 ;i < genotype.size() ;i++ )
     {
         r = rand() % 2;
         if (r)
            genotype[i] = 1;
         else
            genotype[i] = 0;
     }
     point_sub_pop(0)->add_genotype(genotype, N);
     
     if (track_genealogy == 2)
     {
         cout << "initialization generation..." << endl;
         submit_pop_genealogy();
         genealogy.add_generation(max_fitness());
     }
     set_global_generation(0);
     return 0;
}
/**
  * @brief Function to enable tracking global genealogy of the multiple populations.
  * This is overloaded version of the haploid_highd::track_locus_genealogy().
  * If only genealogy of a local sub-population is needed, haploid_highd::track_locus_genealogy() should be used instead.
  * @param vector <int> loci loci to track
  */
int multi_population::track_locus_genealogy(vector <int > loci)
{
    if(HP_VERBOSE){cerr <<"MULTI_POPULATION::track_locus_genealogy(vector <int> loci)... number of loci="<<loci.size();}
    genealogy.reset();
    for (unsigned int i = 0; i < loci.size(); i ++)
    {
        vector <node_t> temp_generation;
        genealogy.newGenerations.push_back(temp_generation);
        genealogy.track_locus(loci[i]);
    }
    genealogy.extend_storage(10 + N());
    if (HP_VERBOSE){cerr<<"done\n";}

    //set_global_generation(-1);
    for (int location_No = 0; location_No < number_of_locations; location_No ++)

    {
        point_sub_pop(location_No)->track_locus_genealogy(loci, 2);
    }
    track_genealogy = 2;
    return 0;
}

/**
  * @brief submits all local genealogies to the global storage
  */
int multi_population::submit_pop_genealogy()
{
    genealogy.newGenerations.clear();
    vector <node_t> temp_generation;
    for (int i = 0; i < genealogy.trees.size(); i ++)
    {
        genealogy.newGenerations.push_back(temp_generation);

    }
    for (int i = 0; i < number_of_locations; i ++)
    {

        for (unsigned int loc_No = 0; loc_No < genealogy.trees.size(); loc_No ++)
        {
            for (unsigned int node_No = 0; node_No <= sub_population[i]->last_clone; node_No ++)
            {
                if (sub_population[i]->genealogy.newGenerations[loc_No][node_No].clone_size != 0)
                {
                    genealogy.newGenerations[loc_No].push_back(sub_population[i]->genealogy.newGenerations[loc_No][node_No]);
                    genealogy.newGenerations[loc_No].back().own_key.location = i;
                }
            }
        }
    }
    return 0;
}
/**
  * @brief returns maximal overall fitness of the population
  @return maximal fitness of the population
  */
double multi_population::max_fitness()
{
    double return_value = 0;
    for (int i = 0; i < number_of_locations; i ++)
    {
        if (return_value < sub_population[i]->get_max_fitness())
            return_value = sub_population[i]->get_max_fitness();

    }
    return return_value;
}


/**
 * @brief Migrate random individual from specified source to a random destination location
 * @param source number of source location
 */
int multi_population::migrate(int source)
{
    int actual_no_migrants = determine_number_of_migrants(source);
    if (actual_no_migrants > 0)
    {
        //cout << "Migrants from location "<< i << ":  "<< actual_no_migrants << endl;
        sub_population[source]->produce_random_sample(actual_no_migrants);
        for (int migrant = 0; migrant < actual_no_migrants; migrant ++)
        {
            int migrant_clone_No = sub_population[source]->random_clone();
            int migrant_clone_destination = determine_migration_destination();
            if (source != migrant_clone_destination){

                transfer_clone(source, migrant_clone_destination, migrant_clone_No);
                number_of_migration_events ++;
            }
        }
        sub_population[source]->random_sample.clear();
    }
    return 0;
}

/**
 * @brief Migrate random individual from a random source to a random destination location
 * @param source number of source location
 */
int multi_population::migrate()
{
    for (int i = 0; i < number_of_locations; i ++)
    {
        int actual_no_migrants = determine_number_of_migrants(i);
        if (actual_no_migrants > 0)
        {
            sub_population[i]->produce_random_sample(actual_no_migrants);
            for (int migrant = 0; migrant < actual_no_migrants; migrant ++)
            {
                int migrant_clone_No = sub_population[i]->random_clone();
                int migrant_clone_destination = determine_migration_destination();
                cout << "Source " << i << " Dest: " << migrant_clone_destination << endl;
                if (i != migrant_clone_destination)
                {
                    transfer_clone(i, migrant_clone_destination, migrant_clone_No);
                    number_of_migration_events ++;
                }
            }
            sub_population[i]->random_sample.clear();
        }
    }
    return 0;
}


/**
  @brief Determine number of migratinf individuals
  The number is derived form the migration rate and the size of a local sub-population
  @param sub_pop_no source location to migrate from
  */
int multi_population::determine_number_of_migrants(int sub_pop_No)
{
<<<<<<< Updated upstream
    int number_of_migrants = gsl_ran_poisson(sub_population[sub_pop_No]->evo_generator, sub_population[sub_pop_No]->N()*migration_rate);
    number_of_migrants = number_of_migrants > MAX_MIGRATION_RATE * point_sub_pop(sub_pop_No)->N() ?
                MAX_MIGRATION_RATE * point_sub_pop(sub_pop_No)->N() : number_of_migrants;
=======
    int number_of_migrants;
    number_of_migrants = gsl_ran_poisson(sub_population[sub_pop_No]->evo_generator, sub_population[sub_pop_No]->N()*migration_rate);
>>>>>>> Stashed changes
    return number_of_migrants;
}

/**
  @brief Finds random location to migrate to
  */
int multi_population::determine_migration_destination()
{
  int destination;
  if (number_of_locations == 1) destination = 0;
  else destination = gsl_rng_uniform_int(sub_population[0]->evo_generator, number_of_locations);
  return destination;

}

/**
  @brief Pick individuum and rasfers it into another location
  This function takes an individum with specified index at the specified location, finds available index in the destination location and transfers clone there. Genealogy is up-dated if needed.
  @param sub_pop_source source location
  @param sub_pop_destination destination location
  @param source index of the clone from which the individuum will be taken
  @return 0
  */
int multi_population::transfer_clone(int sub_pop_source, int sub_pop_destination, int source)
{
    if (sub_population[sub_pop_source]->last_clone < source)
        return -1;
    if (sub_population[sub_pop_source]->population[source].clone_size < 1)
        return -1;


    if (sub_population[sub_pop_destination]->available_clones.size() == 0)
        sub_population[sub_pop_destination]->provide_at_least(1);
    int new_gt = sub_population[sub_pop_destination]->available_clones.back();
    sub_population[sub_pop_destination]->available_clones.pop_back();

    sub_population[sub_pop_destination]->population[new_gt].genotype = sub_population[sub_pop_source]->population[source].genotype;
    sub_population[sub_pop_destination]->population[new_gt].clone_size = 1;
    //Update fitness of the new clone
    sub_population[sub_pop_destination]->calc_individual_traits(sub_population[sub_pop_destination]->population[new_gt]);
    sub_population[sub_pop_destination]->calc_individual_fitness_from_traits(sub_population[sub_pop_destination]->population[new_gt]);
    sub_population[sub_pop_destination]->check_individual_maximal_fitness(sub_population[sub_pop_destination]->population[new_gt]);
    sub_population[sub_pop_destination]->update_traits();
    sub_population[sub_pop_destination]->update_fitness();
    //location[destination]->population_size ++;
    sub_population[sub_pop_destination]->last_clone = (new_gt < sub_population[sub_pop_destination]->last_clone)?sub_population[sub_pop_destination]->last_clone:new_gt;
    //Update general sub-population parameters
    sub_population[sub_pop_source]->population_size --;
    sub_population[sub_pop_source]->population[source].clone_size -- ;
    if (sub_population[sub_pop_source]->population[source].clone_size < 1)
    {
        sub_population[sub_pop_source]->population[source].clone_size = 0;
        sub_population[sub_pop_source]->available_clones.push_back(source);
        sub_population[sub_pop_source]->number_of_clones -- ;
    }
    sub_population[sub_pop_destination]->number_of_clones ++;
    sub_population[sub_pop_destination]->population_size ++;
    // Update parameters for tree
    if (track_genealogy == 2) {
        for (unsigned int genlocus=0; genlocus < genealogy.loci.size(); genlocus++)
        {
           add_migrating_clone_to_genealogy
                    (   genlocus,
                        sub_pop_source,
                        sub_pop_destination,
                        new_gt,
                        source,
                        sub_population[sub_pop_source]->genealogy.newGenerations[genlocus][source].crossover[0],
                       sub_population[sub_pop_source]->genealogy.newGenerations[genlocus][source].crossover[1], 1, 1
                    );
            sub_population[sub_pop_source]->genealogy.newGenerations[genlocus][source].clone_size --;
        }
    }
    return 0;

}


void multi_population::set_global_generation(int new_generation)
{
    generation = new_generation;
    for (int i = 0; i < number_of_locations; i ++)
    {
            sub_population[i]->generation = new_generation;
    }
    return;
}

/**
  @brief Creating record for the migrating clone in the global genealogy
  @param locusIndex index of the genealogy tree to
  @param old_location location of the parent
  @param new_location new location of the migrating individuum
  @param dest index of the new clone in the population vector
  @param parent index of the parent clone in population vector of the source location
  @param left recombination parameter
  @param right recombination parameter
  @param cs recombination parameter
  @param n number of offsprings
  */

void multi_population::add_migrating_clone_to_genealogy(int locusIndex, int old_location, int new_location,  int dest, int parent, int left, int right, int cs, int n)
{

    if (1) {
        cerr <<"multi_population::add_migrating_clone_to_genealogy(): dest:  "<<dest<<" parent: "<<parent<<endl;
        tree_key_t temp;
        temp = sub_population[old_location]->genealogy.newGenerations[locusIndex][parent].parent_node;
        /*temp.age = generation - 1;
        temp.index = sub_population[old_location].newGenerations[locusIndex][parent].parent_node.index;
        temp.location = sub_population[old_location].newGenerations[locusIndex][parent].parent_node.location;*/
        if (genealogy.trees[locusIndex].check_node(temp)){
            cerr <<"multi_population::add_clone_to_genealogy(): parent node ok"<<endl;
        }else{
            cerr <<"multi_population::add_clone_to_genealogy(): parent node DOES NOT EXIST!"<<endl;
        }
    }
    sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].parent_node = sub_population[old_location]->genealogy.newGenerations[locusIndex][parent].parent_node;

    sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].own_key.index = dest;
    sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].own_key.age = generation;
    sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].own_key.location = new_location;

    sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].fitness = sub_population[new_location]->population[dest].fitness;
    sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].number_of_offspring = n;
    sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].clone_size = cs;
    sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].crossover[0] = left;
    sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].crossover[1] = right;

    if (sub_population[new_location]->number_of_traits != sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].traits.size())
    {
        sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].traits.clear();
        for (int traitNo = 0; traitNo < sub_population[new_location]->number_of_traits; traitNo ++)
        {
            sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].traits.push_back(sub_population[new_location]->population[dest].trait[traitNo]);
        }
    }else
    {
        for (int traitNo = 0; traitNo < sub_population[new_location]->number_of_traits; traitNo ++)
        {
            sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].traits[traitNo] = sub_population[new_location]->population[dest].trait[traitNo];

        }
    }

     sub_population[new_location]->genealogy.newGenerations[locusIndex][dest].genotype = sub_population[new_location]->population[dest].genotype;

<<<<<<< Updated upstream


=======
>>>>>>> Stashed changes
    if (HP_VERBOSE) {
        cerr <<"multi_population::add_clone_to_genealogy(): done"<<endl;
    }
    return;
}

/**
  * @brief mutate a sub-population in a specified location
  * This function is overloaded version of haploid_highd::mutate()
  * @param location number
  @return 0
  */
 int multi_population::mutate(int location)
{

    if (HP_VERBOSE)	cerr <<"multi_population::mutate() ..."<<endl;

    vector <int> mutations;
    int tmp_individual=0, nmut=0;
    size_t mutant;
    sub_population[location]->allele_frequencies_up_to_date = false;
    int actual_n_o_mutations, actual_n_o_mutants;
    if (sub_population[location]->mutation_rate > HP_NOTHING and not sub_population[location]->all_polymorphic) {
        //determine the number of individuals that are hit by at least one mutation
        actual_n_o_mutants = gsl_ran_poisson(sub_population[location]->evo_generator, (1.0-exp(-sub_population[location]->mutation_rate*sub_population[location]->number_of_loci))*sub_population[location]->population_size);
        sub_population[location]->produce_random_sample(1 + min(actual_n_o_mutants, sub_population[location]->population_size));

        //make sure enough empty clones are available to accomodate the new mutants
        sub_population[location]->provide_at_least(min(actual_n_o_mutants, sub_population[location]->population_size));

        //loop over the mutant individuals and introduce the mutations
        for (int individual = 0; individual != actual_n_o_mutants; individual++) {
            //determine the target clone
            mutant = sub_population[location]->random_clone();
            //determine the number of mutation it suffers. this should be Poisson conditional on having at least one
            //in practice the solution is fine but it is somewhat inaccurate for multiple mutations
            actual_n_o_mutations = gsl_ran_poisson(sub_population[location]->evo_generator, sub_population[location]->number_of_loci * sub_population[location]->mutation_rate)+1;
            //introduce the mutations, not that flip_single_locus returns the number of new mutant, which is fed back into
            //flip_single_locus to introduce the next mutation
            for (int i = 0; i != actual_n_o_mutations; i++)
                mutant=flip_single_locus(location, mutant, gsl_rng_uniform_int(sub_population[location]->evo_generator, sub_population[location]->number_of_loci));
        }
    } else if(sub_population[location]->all_polymorphic) {
        if(HP_VERBOSE) cerr <<"haploid_highd::mutate(): keeping all loci polymorphic"<<endl;
        sub_population[location]->calc_allele_freqs(); //calculate the allele frequencies
        nmut=0;
        for (int locus=0; locus<sub_population[location]->L(); locus++){	//loop over all loci
            if (fabs(2*sub_population[location]->allele_frequencies[locus]-1)>1-HP_NOTHING){	//spot fixed loci
                if ((sub_population[location]->ancestral_state[locus]==0 and (2*sub_population[location]->allele_frequencies[locus]-1)<0) or
                    (sub_population[location]->ancestral_state[locus]==1 and (2*sub_population[location]->allele_frequencies[locus]-1)>0))
                {	//if they are in the ancestral state
                    tmp_individual = flip_single_locus(location, locus);		//introduce new allele
                    sub_population[location]->polymorphism[locus].birth = sub_population[location]->get_generation();
                    sub_population[location]->polymorphism[locus].fitness = sub_population[location]->population[tmp_individual].fitness-sub_population[location]->fitness_stat.mean;
                    sub_population[location]->polymorphism[locus].fitness_variance = sub_population[location]->fitness_stat.variance;
                    nmut++;
                }else{	//if locus is in derived state, flip coefficient of trait zero
                    sub_population[location]->trait[0].set_additive_coefficient(-sub_population[location]->trait[0].get_additive_coefficient(locus),locus,locus);
                    sub_population[location]->fixed_mutations.push_back(sub_population[location]->polymorphism[locus]);
                    sub_population[location]->fixed_mutations.back().sweep_time = sub_population[location]->get_generation() -sub_population[location]->fixed_mutations.back().birth;
                    tmp_individual=flip_single_locus(location, locus);
                    sub_population[location]->ancestral_state[locus]= (sub_population[location]->ancestral_state[locus]==0)?1:0;
                    sub_population[location]->polymorphism[locus].birth = sub_population[location]->get_generation();
                    sub_population[location]->polymorphism[locus].effect = (2*sub_population[location]->ancestral_state[locus]-1)*sub_population[location]->trait[0].get_additive_coefficient(locus);
                    sub_population[location]->polymorphism[locus].fitness = sub_population[location]->population[tmp_individual].fitness;
                    sub_population[location]->polymorphism[locus].fitness_variance = sub_population[location]->fitness_stat.variance;
                    nmut++;
                }
            }
        }
        sub_population[location]->number_of_mutations.push_back(nmut);
        sub_population[location]->calc_stat();

    } else if(HP_VERBOSE) cerr <<"haploid_highd::mutate(): mutation rate is zero."<<endl;

    if (HP_VERBOSE)	cerr <<"done."<<endl;;
    return 0;
}


/**
  @brief A single mutation event of a specified individuum in a specified location\
  @param location
  @param clonenum
  @param locus
  @return index of the new clone
  */
unsigned int multi_population::flip_single_locus(int location, unsigned int clonenum, int locus) {
    // produce new genotype
    int new_clone = sub_population[location]->available_clones.back();
    sub_population[location]->available_clones.pop_back();
    sub_population[location]->allele_frequencies_up_to_date = false;

    //copy old genotype
    sub_population[location]->population[new_clone].genotype = sub_population[location]->population[clonenum].genotype;
    // new clone size == 1, old clone reduced by 1
    sub_population[location]->population[new_clone].clone_size = 1;
    sub_population[location]->population[clonenum].clone_size--;
    // flip the locus in new clone
    sub_population[location]->population[new_clone].genotype.flip(locus);
    // calculate traits and fitness
    vector<int> diff(1, locus);
    for (int t = 0; t < sub_population[location]->number_of_traits; t++){
        sub_population[location]->population[new_clone].trait[t] = sub_population[location]->population[clonenum].trait[t] + sub_population[location]->get_trait_difference(sub_population[location]->population[new_clone], sub_population[location]->population[clonenum], diff, t);
    }
    sub_population[location]->calc_individual_fitness_from_traits(sub_population[location]->population[new_clone]);
    sub_population[location]->check_individual_maximal_fitness(sub_population[location]->population[new_clone]);

    //update the last clones that is to be tracked
    sub_population[location]->last_clone = (new_clone<sub_population[location]->last_clone)?sub_population[location]->last_clone:new_clone;

    // add clone to current population
    if (sub_population[location]->population[clonenum].clone_size == 0)
        sub_population[location]->available_clones.push_back(clonenum);
    else
        sub_population[location]->number_of_clones++;

    if (track_genealogy == 2) {
        for (unsigned int genlocus=0; genlocus<genealogy.loci.size(); genlocus++) {
            add_migrating_clone_to_genealogy(
                    genlocus, location, location, new_clone, clonenum,
                    sub_population[location]->genealogy.newGenerations[genlocus][clonenum].crossover[0],
                    sub_population[location]->genealogy.newGenerations[genlocus][clonenum].crossover[1], 1, 1
                    );
            sub_population[location]->genealogy.newGenerations[genlocus][clonenum].clone_size--;
        }
    }


    if (HP_VERBOSE >= 2) cerr <<"subpop::flip_single_spin(): mutated individual in clone "<<clonenum<<" at locus "<<locus<<endl;
    return new_clone;
}

/**
  @brief A single mutation event of a randomly chosen individuum in a specified location
  @param location
  @param locus
  */
unsigned int multi_population::flip_single_locus(int location, int locus) {
    if (sub_population[location]->available_clones.size() == 0)
        sub_population[location]->provide_at_least(1);
    return flip_single_locus(location, sub_population[location]->random_clone(), locus);
}

/**
  @brief perform evolve cycle on the population
  This function is the overloaded version of the haploid_highd::evolve(int gen)
  @param gen number of generations to evolve
  @return 0 if OK, number of errors occured otherwise;
  */
int multi_population::evolve(int gen)
{
    int err = 0;
    int g = 0;
    while (g < gen && err == 0)
    {
        for (int cur_loc = 0 ; cur_loc < number_of_locations; cur_loc ++)
        {
            if (sub_population[cur_loc]->N() > 0)
            {
                sub_population[cur_loc]->update_traits();
                sub_population[cur_loc]->update_fitness();
                //evolve of the haploid highd (not used)
                //err +=  sub_population[cur_loc].evolve_loc(cur_loc, gen);
                //evolve() of the multi_pop
                err += evolve_local(cur_loc, 1);
            }
        }
        if (N() == 0)
        {
<<<<<<< Updated upstream
            //TODO throw population went extinct exception!!!
            cerr << "Unfortunately, poplation went extinct!" << endl;
            throw HP_EXTINCTERR;
        }
        migrate();
        generation ++;
        g++;
        set_global_generation(generation);

        //add the current generation to the genealogies and prune (i.e. remove parts that do not contribute the present.
        for (int location = 0; location < number_of_locations; location ++)
        {
           if (track_genealogy == 1) sub_population[location]->genealogy.add_generation(fitness_max);
        }
        if (track_genealogy == 2)
        {
            submit_pop_genealogy();
            genealogy.add_generation(max_fitness());
        }
=======
            sub_population[cur_loc]->update_traits();
            sub_population[cur_loc]->update_fitness();
            //evolve of the haploid highd (not used)
            //err +=  sub_population[cur_loc].evolve_loc(cur_loc, gen);
            //evolve() of the multi_pop
            err += evolve_local(cur_loc, gen);

        }
    }
    //migrate();


    //add the current generation to the genealogies and prune (i.e. remove parts that do not contribute the present.
    for (int location = 0; location < number_of_locations; location ++)
    {
       if (track_genealogy == 1) sub_population[location]->genealogy.add_generation(fitness_max);
>>>>>>> Stashed changes
    }
//    if (track_genealogy == 2)
//    {
//        submit_pop_genealogy();
//        genealogy.add_generation(max_fitness());
//    }
    //generation ++;
    //set_global_generation(generation);
    cout << population_size() << endl;
    return err;
}

<<<<<<< Updated upstream

=======
/**
  * @brief perform evolve cycle on a sub-population in a specified location
  * @param location location number
  * @param gen number of generations to evolve
  @return 0 if OK, number of errors occured otherwise;

  */
>>>>>>> Stashed changes
int multi_population::evolve_local(int location, int gen) {
    if (HP_VERBOSE) cerr<<"multi_population::evolve(int gen)...";
    if (sub_population[location]->population_size == 0)
        return 0;
    if (N() == 0)
    {
        cout << "Population went extinct!" << endl;
        return 0;
    }
    int err=0, g=0;
    sub_population[location]->allele_frequencies_up_to_date = false;
    // calculate an effective outcrossing rate to include the case of very rare crossover rates.
    // Since a recombination without crossovers is a waste of time, we scale down outcrossing probability
    // and scale up crossover rate so that at least one crossover is guaranteed to happen.
     if (sub_population[location]->recombination_model == CROSSOVERS)
        sub_population[location]->outcrossing_rate_effective = sub_population[location]->outcrossing_rate * (1 - exp(-sub_population[location]->number_of_loci * sub_population[location]->crossover_rate));
    else
        sub_population[location]->outcrossing_rate_effective = sub_population[location]->outcrossing_rate;

    // evolve cycle
    while((err == 0) && (g < gen)) {
        if (1) cerr << "generation " << generation << endl;
        sub_population[location]->random_sample.clear();			//discard the old random sample
        if(err==0) err=sub_population[location]->select_gametes();	//select a new set of gametes (partitioned into sex and asex)
        else if(HP_VERBOSE) cerr<<"Error in select_gametes()"<<endl;
        sort(sub_population[location]->available_clones.begin(), sub_population[location]->available_clones.end(), std::greater<int>()); //sort clones in order to use the first ones again and again
        //if(err==0) err=sub_population[location]->add_recombinants();	//do the recombination between pairs of sex gametes
        //else if(HP_VERBOSE) cerr<<"Error in recombine()"<<endl;
        if(err==0 && sub_population[location]->N() > 0) err=mutate(location);		//mutation step
        else if(HP_VERBOSE) cerr<<"Error in mutate()"<<endl;
        sub_population[location]->random_sample.clear();			//discard the old random sample
        g++;
        //sub_population[location]->generation++;
        //FIXME increment generation?
        //generation ++;
    }
    if (1) {
        if(err==0) cerr<<"done."<<endl;
        else cerr<<"error "<<err<<"."<<endl;
    }
    return err;
}



void multi_population::set_mutation_rate(double rate){
    for (int i = 0; i < number_of_locations; i ++)
    {

        point_sub_pop(i)->set_mutation_rate(rate);
    }
}

void multi_population::set_carrying_capacity(int capacity){
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->carrying_capacity = capacity;
    }
}

void multi_population::set_outcrossing_rate(double o_rate){
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->outcrossing_rate = o_rate;
    }
}

void multi_population::set_crossover_rate(double c_rate){
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->crossover_rate = c_rate;
    }
}

void multi_population::set_recombination_model(int r_model){
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->recombination_model = FREE_RECOMBINATION;
    }
}

void multi_population::set_trait_coefficient(double coefficient, vector<int> loci, int trait_no){
    for (int i = 0; i < number_of_locations; i ++)
    {
        point_sub_pop(i)->add_trait_coefficient(coefficient, loci, trait_no);
        point_sub_pop(i)->update_traits();
        point_sub_pop(i)->update_fitness();
    }
}

void multi_population::set_trait_weights(double* weights){
    for (int i = 0; i < number_of_locations; i ++){
        point_sub_pop(i)->set_trait_weights(weights);
        point_sub_pop(i)->update_traits();
        point_sub_pop(i)->update_fitness();
    }

}

void multi_population::add_random_genotype(int N_in){
    boost::dynamic_bitset<> genotype(L);
    bool r = 0;
    for( int i = 0 ;i < genotype.size() ;i ++ )
    {
        r = rand() % 2;
        if (r)
           genotype[i] = 1;
        else
           genotype[i] = 0;
    }
    point_sub_pop(0)->add_genotype(genotype, N_in); //FIXME we can only set genotype to location 0, otherwise genealogy gets crazy.
    set_global_generation(0);
    if (track_genealogy == 2){
        submit_pop_genealogy();
        genealogy.add_generation(max_fitness());
    }
}
