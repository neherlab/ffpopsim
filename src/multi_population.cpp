#include "multi_population.h"

//#define SELECT 0.01
class hh_source : public haploid_highd{

    //Center of the distribution in polar coordinates:
    //double phi_0 = 0;
    //double rho_0 = 0;
    //VARIANCE
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





multi_population::multi_population(int new_locations, int L_in, int n_o_traits, int rng_seed)
{

    //sub_population = new haploid_highd[new_locations]; //.push_back(haploid_highd());


    //Location number check:
    if (new_locations < 1)
    {
        if (HP_VERBOSE)
            cerr << "too few locations provided for the required type of simulations!"<< endl;
        throw(HP_BADARG);
    }

    try{
        environmental_hh* loc;
        for (int i = 0; i < new_locations; i ++)
        {
            loc = new environmental_hh();
            loc->set_up(L_in, rng_seed, n_o_traits);
            sub_population.push_back(loc);
        }
    }catch (int err){
            throw err;
    }
    number_of_traits = n_o_traits;
    MAX_MIGRATION_RATE = 0.5;
    number_of_loci = L_in;
    generation = 0;
    number_of_locations = new_locations;
    mutation_rates.resize(number_of_locations);
    carrying_capacities.resize(number_of_locations);

    vector <double> temp;
    for (int i = 0; i < number_of_locations; i ++){
        temp.push_back(0);
    }

    for (int i = 0; i < number_of_locations; i ++){
        migration_rates.push_back(temp);
    }


    track_genealogy = 0; //No genealogy by default
    critical_migration_rate = 0.0; //No migration by default
    number_of_migration_events = 0;
    //population_sizes.resize(new_locations);
    for (int i = 0; i < new_locations; i ++ ){
        //population_sizes[i] = 0;
    }


}

multi_population::~multi_population()
{
    for (int i = 0; i < number_of_locations; i ++)
     {
        //population_sizes.clear();
        sub_population.clear();
        //delete sub_population[i];
     }
}

int multi_population::N(int i){
    if (i == -1){
        int size = 0;
        for (int loc = 0; loc < number_of_locations; loc ++){
            size = size + point_sub_pop(0 )->N();
        }
        return size;
    }else{
        if (i < 0 || i > number_of_locations){
            cout << "No location found!!" << endl;
            throw (HP_BADARG);
        }
        return point_sub_pop(i)->N();
    }
}

void multi_population::reset()
{
    number_of_locations = 0;
    sub_population.clear();
    //population_sizes.clear();

}

void multi_population::track_locus_genealogy(vector <int > loci)
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
    return;
}

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

double multi_population::max_fitness()
{
    double return_value = 0;
    for (int i = 0; i < number_of_locations; i ++)
    {
        if (return_value > sub_population[i]->get_max_fitness())
            return_value = sub_population[i]->get_max_fitness();
    }
    return return_value;
}


/**
 @brief multi_population::migrate
 * @return
 */



void multi_population::set_critical_migration_rate(double new_rate){
    if (new_rate < 0 || new_rate > 1){
        cout << "Do you really want to set this value?" << endl;
        throw(HP_BADARG);
    }
    critical_migration_rate = new_rate;
}

void multi_population::set_migration_rates(vector<double> new_migration_rates, int source){
    if (new_migration_rates.size() != number_of_locations ||
            source < 0 ||
            source >= number_of_locations){
        cout << "Bad argument! Either array size mismatch, or location does not exist" << endl;
        throw(HP_BADARG);
    }
    double total_m_r = 0;
    for (int i = 0; i < number_of_locations; i ++){
        total_m_r += new_migration_rates[i];
    }
    if (total_m_r > critical_migration_rate){
        cout << "Migration values exceed critical migration rate of the population!!" << endl;
        throw(HP_BADARG);
    }
    for (int i = 0; i < number_of_locations; i ++){
        migration_rates[source][i] = new_migration_rates[i];
    }
    return;
}

int multi_population::set_migration_rates(vector< vector<double> > new_migration_rates){
    if (new_migration_rates.size() != number_of_locations){
        cout << "" << endl;
        throw(HP_BADARG);
    }
    for (int i = 0; i < number_of_locations; i ++){
        if (new_migration_rates[i].size() != number_of_locations){
            cout << "The size of the migration matrix must be (locations X locations) !!" << endl;
            throw(HP_BADARG);
        }else{
            set_migration_rates(new_migration_rates[i], i);
        }
    }

    return 0;
}

void multi_population::migrate(){
    vector <int> pop_sizes;
    for (int i = 0; i < number_of_locations; i ++){
        pop_sizes.push_back(point_sub_pop(i)->N());
    }
    for (int source = 0; source < number_of_locations; source ++){
        for (int destination = 0; destination < number_of_locations; destination ++){
            double m_r = migration_rates[source][destination];//particular rate of migration between the locations
            int migrants_num = pop_sizes[source] * m_r;
            if (migrants_num > 0 && source != destination){
                for (int clone_num = 0; clone_num < migrants_num; clone_num ++){
                    int migrant_clone_No = sub_population[source]->random_clone();
                    transfer_clone(source, destination, migrant_clone_No);
                    number_of_migration_events ++;
                }
            }
        }
    }
    return;
}

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

void multi_population::add_migrating_clone_to_genealogy(int locusIndex, int old_location, int new_location,  int dest, int parent, int left, int right, int cs, int n)
{

    if (HP_VERBOSE) {
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





    if (HP_VERBOSE) {
        cerr <<"multi_population::add_clone_to_genealogy(): done"<<endl;
    }
    return;
}

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

unsigned int multi_population::flip_single_locus(int location, int locus) {
    if (sub_population[location]->available_clones.size() == 0)
        sub_population[location]->provide_at_least(1);
    return flip_single_locus(location, sub_population[location]->random_clone(), locus);
}

int multi_population::evolve(int gen)
{
    cout << N() << endl;
    int err = 0;
    int g = 0;
    while (g < gen && err == 0)
    {
        for (int cur_loc = 0 ; cur_loc < number_of_locations; cur_loc ++){
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
        if (N() == 0){
            cerr << "Unfortunately, poplation went extinct!" << endl;
            throw HP_EXTINCTERR;
        }
        migrate();
        generation ++;
        g++;
        set_global_generation(generation);
        update_population_size();
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
    }
}

int multi_population::evolve_local(int location, int gen) {
    if (HP_VERBOSE) cerr<<"multi_population::evolve(int gen)...";
    if (sub_population[location]->population_size == 0)
        return 0;
    if (N() == 0)
    {
        cerr << "Population went extinct!" << endl;
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
        if (HP_VERBOSE) cerr << "generation " << generation << endl;
        sub_population[location]->random_sample.clear();			//discard the old random sample
        if(err==0) err=sub_population[location]->select_gametes();	//select a new set of gametes (partitioned into sex and asex)
        else if(HP_VERBOSE) cerr<<"Error in select_gametes()"<<endl;
        sort(sub_population[location]->available_clones.begin(), sub_population[location]->available_clones.end(), std::greater<int>()); //sort clones in order to use the first ones again and again
        if(err==0) err=sub_population[location]->add_recombinants();	//do the recombination between pairs of sex gametes
        else if(HP_VERBOSE) cerr<<"Error in recombine()"<<endl;
        if(err==0 && sub_population[location]->N() > 0) err=mutate(location);		//mutation step
        else if(HP_VERBOSE) cerr<<"Error in mutate()"<<endl;
        sub_population[location]->random_sample.clear();			//discard the old random sample
        g++;
        //sub_population[location]->generation++;
        //generation ++;
    }
    cout << N()<< endl;
    if (HP_VERBOSE) {
        if(err==0) cerr<<"done."<<endl;
        else cerr<<"error "<<err<<"."<<endl;
    }
    return err;
}

void multi_population::set_mutation_rates(vector<double> rates){
    if (rates.size() != number_of_locations){
        if (HP_VERBOSE){
            cout << "The input vector size must match thre number of the locations!" << endl;
        }
        throw (HP_BADARG);
    }
    mutation_rates = rates;
    for (int i = 0; i < number_of_locations; i ++)
    {

        point_sub_pop(i)->set_mutation_rate(mutation_rates[i]);
    }

}

double multi_population::get_mutation_rate(int location){
    if (location < 0 || location >= number_of_locations){
        if (HP_VERBOSE){
            cout << "Cannot return mutation rate for the location: "<< location << ". No locaiton found!" << endl;
        }
        throw (HP_BADARG);
    }
    return mutation_rates[location];
}

void multi_population::set_carrying_capacities(vector<int> capacities){
    if (capacities.size() != number_of_locations){
        if (HP_VERBOSE){
            cout << "The input vector size must match thre number of the locations!" << endl;
        }
        throw (HP_BADARG);
    }
    carrying_capacities = capacities;
    for (int i = 0; i < number_of_locations; i ++){
        point_sub_pop(i)->carrying_capacity = carrying_capacities[i];
    }
}

double multi_population::get_carrying_capacity(int location){
    if (location < 0 || location >= number_of_locations){
        if (HP_VERBOSE){
            cout << "Cannot return carrying capacity for the location: "<< location << ". No locaiton found!" << endl;
        }
        throw (HP_BADARG);
    }
    return carrying_capacities[location];
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
    for (int i = 0; i < number_of_locations; i ++){
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

clone_t multi_population::get_random_clone(int location){

    int index = 0;
    int clone_index = 0; //point_sub_pop(index)->random_clone();

    if (location == -1){
        if (number_of_locations > 1){
            //get random location
            vector<int> non_empty_locations;
            for (int i = 0; i < number_of_locations; i ++){
                if (point_sub_pop(i)->N() > 0)
                    non_empty_locations.push_back(i);
            }
            if (non_empty_locations.size() > 1){
                index = gsl_ran_flat(point_sub_pop(0)->evo_generator, 0, non_empty_locations.size());
                clone_index = point_sub_pop(non_empty_locations[index])->random_clone();
                return point_sub_pop(non_empty_locations[index])->population[clone_index];
            }else{
                clone_index = point_sub_pop(non_empty_locations[0])->random_clone();
                return point_sub_pop(non_empty_locations[0])->population[clone_index];
            }
        }else{
            clone_index = point_sub_pop(0)->random_clone();
            return point_sub_pop(0)->population[clone_index];
        }
    }else{
        if (location > number_of_locations){
            cout << "Location with this number does not exist!" << endl;
            throw HP_BADARG;
        }
        if (point_sub_pop(location)->N() < 1){
            cout << "Location spot is empty! returning dummy clone" << endl;


        }
        clone_index = point_sub_pop(location)->random_clone();
        return point_sub_pop(location)->population[clone_index];
    }
}

boost::dynamic_bitset<> multi_population::get_genotype(clone_t clone){
    return clone.genotype;
}


void multi_population::set_random_genotype(int N_in){
    boost::dynamic_bitset<> genotype(number_of_loci);
    bool r = 0;
    for( uint i = 0 ;i < genotype.size() ;i ++ )
    {
        r = rand() % 2;
        if (r)
           genotype[i] = 1;
        else
           genotype[i] = 0;
    }
    for (int i = 0; i < number_of_locations; i ++ ){
        point_sub_pop(i)->allele_frequencies_up_to_date = false;
    // Clear population
    point_sub_pop(i)->population.clear();
    point_sub_pop(i)->available_clones.clear();
    if (track_genealogy) {
        point_sub_pop(i)->genealogy.reset_but_loci();
    }

    point_sub_pop(i)->population_size = 0;
    point_sub_pop(i)->random_sample.clear();

    // Initialize the clones and calculate the population size
    point_sub_pop(i)->number_of_clones = 0;
    point_sub_pop(i)->last_clone = 0;
    point_sub_pop(i)->provide_at_least(N_in);
    }
    point_sub_pop(0)->add_genotype(genotype, N_in); //FIXME we can only set genotype to location 0, otherwise genealogy gets crazy.
    //set_global_generation(0); //FIXME
    // set the carrying capacity if unset
    if(get_carrying_capacity(0) < HP_NOTHING){
        vector <int> capacities;
        for(int i = 0; i < number_of_locations; i ++){
            capacities.push_back(N());
        }
        set_carrying_capacities(capacities);
    }
    // Calculate all statistics to be sure
    set_global_generation(generation++);
    for (int i = 0; i < number_of_locations; i ++){
        point_sub_pop(i)->calc_stat();
    }
    //add the current generation to the genealogies and prune (i.e. remove parts that do not contribute the present)
    if (track_genealogy == 1){genealogy.add_generation(fitness_max);}

    if (track_genealogy == 2){
        submit_pop_genealogy();
        genealogy.add_generation(max_fitness());
    }
    if (HP_VERBOSE) cerr <<"done."<<endl;
    update_population_size();
}

void multi_population::set_environment(double rho, double phi, int location){
    if (location == -1){
        for (int loc = 0; loc < number_of_locations; loc ++){
            set_environment(rho,phi,loc);
        }
    }else{
        if (location < 0 || location > number_of_locations){
            throw HP_BADARG;
        }
        point_sub_pop(location)->setRho0(rho);
        point_sub_pop(location)->setPhi0(phi);
    }
}

void multi_population::update_population_size(){
    for (int i = 0; i < number_of_locations; i ++){
        //population_sizes[i] = point_sub_pop(i)->N();
    }

}
void multi_population::update_fitness(){
    for (int loc =0; loc < number_of_locations; loc ++){
        point_sub_pop(loc)->update_fitness();
    }
}
