void multi_population::set_carrying_capacity(int N) {
    for (int i = 0; i < number_of_locations, i ++)
    {
        pop.point_sub_pop(i)->carrying_capacity = N;
    }
}
void multi_population::set_mutation_rate(double mu) {
    for (int i = 0; i < number_of_locations, i ++)
    {
        pop.point_sub_pop(i)->set_mutation_rate(mu);
    }
}
void multi_population::set_outcrossing_rate(double o_rate) {
    for (int i = 0; i < number_of_locations, i ++)
    {
        pop.point_sub_pop(i)->outcrossing_rate = o_rate;
    }
}
void multi_population::set_carrying_capacity_global(double c_rate) {
    for (int i = 0; i < number_of_locations, i ++)
    {
        pop.point_sub_pop(i)->crossover_rate = c_rate;//1e-2;
    }
}

void multi_population::set_recombination_model(int rec_model)
    for (int i = 0; i < number_of_locations, i ++)
    {
        pop.point_sub_pop(i)->recombination_model = rec_model;
    }
}

void multi_population::set_trait_coefficients(int coefficient, vector <int> loci, int trait_no)
{
    for (int i = 0; i < number_of_locations, i ++)
    {
        pop.point_sub_pop(i)->add_trait_coefficient(0, loci, 1);
    }
}

void multi_population::set_trait_weights(double* new_weights)
{
    for (int i = 0; i < number_of_locations, i ++)
    {
        pop.point_sub_pop(i)->set_trait_weights(new_weights);
    }
}

void multi_population::update_traits()
{
    for (int i = 0; i < number_of_locations, i ++)
    {
        pop.point_sub_pop(i)->update_traits();
    }
}
void multi_population::update_fitness()
{
    for (int i = 0; i < number_of_locations, i ++)
    {
        pop.point_sub_pop(i)->update_fitness();
    }
}


void multi_population::add_random_genotype(int N)
{
    boost::dynamic_bitset<> genotype(number_of_loci);
    for( int i = 0 ;i < genotype.size() ;i++ )
     {
         r = rand() % 2;
         if (r)
            genotype[i] = 1;
         else
            genotype[i] = 0;
     }
     pop.point_sub_pop(0)->add_genotype(genotype, N);
     
     if (track_genealogy == 2)
     {
        pop.submit_pop_genealogy();
        pop.genealogy.add_generation(pop.max_fitness());
     }
}

       

