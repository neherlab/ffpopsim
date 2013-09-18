
/* Include directives */
#include "ffpopsim_highd.h"
#include "multi_population.h"
#include <cstring>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cmath>
#include <boost/foreach.hpp>
//#include <boost/filesystem.hpp>
#define HIGHD_BADARG -1354341
#define NOTHING 1e-10

/* Be verbose? */
#define HIGHD_VERBOSE 0


//#define MUTATION_RATE 1e-5
//#define MIGRATION_RATE 0
#define NUMBER_OF_TAITS 2
#define _L 100
#define GENERATIONS_Eq  150000
#define SHIFT_BEGIN 20000
#define SHIFT_END 90000
//#define SELECTIVE 1
//#define SHIFT_RATE 5000
#define GENERATIONS  10000
#define Av_Num 0
//#define SAMPLE_SIZE 1500



void shift_landscape (multi_population* pop, double delta_phi)
{
    for (int i = 0; i < pop->get_locations(); i ++)
    {
        pop->point_sub_pop(i)->phi_0 = pop->point_sub_pop(i)->phi_0 + delta_phi;
    }
}

void shift_landscape_random (multi_population* pop, int location, double sigma = 3.14/8)
{
    double delta_phi = 0; //gsl_ran_gaussian_ziggurat (pop->point_sub_pop(location)->evo_generator, sigma);
    pop->point_sub_pop(location)->phi_0 = pop->point_sub_pop(location)->phi_0 + delta_phi;
    return;
}

void sample_genotypes(multi_population* pop, string genotype_file_name, int time)
{

    ofstream genotype_file;
    //stringstream ss;
    //ss << pop->generation;
    //genotype_file_name.append(ss.str());
    const char * genotype_fn = genotype_file_name.c_str();
    genotype_file.open(genotype_fn, std::ios_base::app);

    for (int loc_num = 0; loc_num < pop->get_locations(); loc_num ++)
    {
        if (pop->point_sub_pop(loc_num)->N() == 0)
        {
            return;
        }
        int sample_size = (pop->point_sub_pop(loc_num)->N() > 100) ? 100 : pop->point_sub_pop(loc_num)->N();

        //produce random sample

        pop->point_sub_pop(loc_num)->random_sample.clear();
        pop->point_sub_pop(loc_num)->produce_random_sample(sample_size);

        double mean_genotype[_L];
        for (int i = 0; i < _L; i ++ )
        {
            mean_genotype[i] = 0;
        }

        for (int sample = 0; sample < sample_size; sample ++)
        {
            for (int pos = 0; pos < _L; pos ++)
            {
                if (pop->point_sub_pop(loc_num)->get_genotype_string(pop->point_sub_pop(loc_num)->random_sample[sample])[pos] == '1')
                {
                    mean_genotype[pos] += 1.0 / sample_size;
                }
            }
        }

        genotype_file << loc_num << "    "<< time << "    ";
        for (int i = 0; i < _L; i ++ )
        {
            genotype_file  << mean_genotype[i] << " ";
        }
        genotype_file << endl;
    }
    genotype_file.close();

}

void sample(multi_population* pop, string traits_file_name)
{

    ofstream traits_file;
    const char * traits_fn = traits_file_name.c_str();
    traits_file.open(traits_fn, std::ios_base::app);

    for (int location = 0 ; location < pop->get_locations(); location ++ )
    {

        stat_t fitness_stat1;
        fitness_stat1 = pop->point_sub_pop(location)->get_fitness_statistics();

        stat_t t_stat1;
        t_stat1 = pop->point_sub_pop(location)->get_trait_statistics(0);
        stat_t t_stat2;
        t_stat2 = pop->point_sub_pop(location)->get_trait_statistics(1);
        stat_t t_stat3;
        t_stat3 = pop->point_sub_pop(location)->get_trait_statistics(2);

        //traits_file << location <<  "   " << pop->generation << "   "<< t_stat1.mean << "  "<<   t_stat1.variance << "     "<<  t_stat2.mean << "  "<< t_stat2.variance << "  "<< t_stat3.mean <<   "  " << t_stat3.variance <<   "  " << fitness_stat1.mean << "  "<< fitness_stat1.variance << "  "<< pop->point_sub_pop(location)->phi_0<<  endl;

    }

    traits_file.close();
}



int main(int argc, char  *argv[]){
//        // ############################################
//        //
//        // Read the input parameters of the simulation
//        //
//        // ############################################
//        int SELECTIVE = atof(argv[1]); //MIGRATION_RATE = atof(argv[1]);
//        double POPULATION_SIZE = atoi(argv[2]);
//        const char * fileDir =  argv[3];
//        int locations = atoi(argv[4]);
//        double MUTATION_RATE = atof(argv[5]);
//        double MIGRATION_RATE = atof(argv[6]);
//        // int _L = atoi(argv[5]);

//        // ############################################

//        // ############################################
//        //
//        // Initialize the essential parameters for the simulation and load the population
//        //
//        // ############################################

//        int N = 20, N_pop = POPULATION_SIZE;
          multi_population pop(10, 100);
          //pop.set_migration_rate(0.1);

          //pop.set_carrying_capacity(100);
//        pop.set_mutation_rate(0.001);
//        pop.set_crossover_rate(0.0);
//        pop.set_outcrossing_rate(0.0);
//        pop.set_recombination_model(0);

          vector <int> gen_loci;
          gen_loci.push_back(6);
          pop.track_locus_genealogy(gen_loci);

          pop.set_random_genotype(20);
        // ############################################


//        // ############################################

//        //#############################################
//        //#
//        //# Set genotype -> trait -> fitness maps
//        //#
//        //#############################################
//        vector <int> loci;
//        for(int j = 0; j < _L ; j++){
//            loci.push_back(j);
//        }
//        for (int trait_no = 0; trait_no < NUMBER_OF_TAITS; trait_no ++){
//            pop.set_trait_coefficient(0.2, loci, trait_no);
//        }
//        /*Set trait weights*/
//        double * weights = new double[NUMBER_OF_TAITS];
//        for (int weight_No = 0; weight_No < NUMBER_OF_TAITS; weight_No ++)
//        {
//            weights[weight_No] = SELECTIVE;
//        }
//        pop.set_trait_weights(weights);
//        delete [] weights;




//        //#############################################
//        //#
//        //# Set the population and track the genealogy for the population
//        //#
//        //#############################################

////        boost::dynamic_bitset<> genotype(_L);
////            bool r = 0;
////            for( int i = 0 ;i < genotype.size() ;i++ )
////            {
////                r = rand() % 2;
////                if (r)
////                   genotype[i] = 1;
////                else
////                   genotype[i] = 0;
////            }
////            pop.point_sub_pop(0)->add_genotype(genotype, 1000);

////        //pop.point_sub_pop(0)->set_wildtype(100);


////        pop.set_global_generation(0);
////        pop.submit_pop_genealogy();
////        pop.genealogy.add_generation(pop.max_fitness());



//        stat_t fitstat;
//        gsl_histogram *SFS = gsl_histogram_alloc(20);
//        gsl_histogram_set_ranges_uniform(SFS,0,1);
//        // ############################################

        //int generation = 0;
        for (int i = 0; i < 1000; i ++)
        {
            //pop.migrate(0);
            pop.evolve(2);
            //pop.migrate(0);


            if (i % 1 == 0) cout << "generation: " << pop.get_generation() << "   " << pop.N() << "    "<<pop.number_of_migration_events<<  endl;
        }
return 0;
}

