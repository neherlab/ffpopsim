
/**
 * @file test_genealogies.cpp
 * @brief Tests for the genealogies in the high-dimensional simulation library.
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version
 * @date 2012-04-20
 */

/* Include directives */
#include "ffpopsim_highd.h"
#include "multi_population.h"
#include <cstring>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cmath>
#include <boost/foreach.hpp>

#define HIGHD_BADARG -1354341
#define NOTHING 1e-10

/* Be verbose? */

#define HIGHD_VERBOSE 0


#define MUTATION_RATE 1e-4
#define NUMBER_OF_TAITS 2
//#define SAMPLE_SIZE 1500
#define _L 100



double get_inter_pop_divergence (multi_population pop, int loc_1, int loc_2, int sample_size)
{
    double result = 0 ;
    if ( pop.get_locations() < loc_1 || pop.get_locations() < loc_2 )
        return 0;
    if ( pop.point_sub_pop(loc_1)->N() < 1 || pop.point_sub_pop(loc_2)->N() < 1 )
        return 0;
    int actual_sample_size = 0;
    if ( pop.point_sub_pop(loc_1)->N() < sample_size || pop.point_sub_pop(loc_2)->N() < sample_size )
        actual_sample_size = min (pop.point_sub_pop(loc_1)->N(), pop.point_sub_pop(loc_2)->N());
    else actual_sample_size = sample_size;

    pop.point_sub_pop(loc_1)->random_sample.clear();
    pop.point_sub_pop(loc_2)->random_sample.clear();
    pop.point_sub_pop(loc_1)->produce_random_sample(actual_sample_size);
    pop.point_sub_pop(loc_2)->produce_random_sample(actual_sample_size);

    for (int i = 0; i < actual_sample_size; i ++)
    {
        string str1 = pop.point_sub_pop(loc_1)->get_genotype_string(pop.point_sub_pop(loc_1)->random_sample[i]);
        string str2 = pop.point_sub_pop(loc_2)->get_genotype_string(pop.point_sub_pop(loc_2)->random_sample[i]);
        int gen_size = str1.size();

        for (int pos = 0; pos < gen_size; pos ++)

        {
            if (str1[pos] != str2[pos])

            {

                result += 1;

            }


        }


    }
    return result / actual_sample_size;
}


double get_segregating_sites_num (multi_population pop, int loc_1, int loc_2, int sample_size)
{
    double result = 0 ;
    if ( pop.get_locations() < loc_1 || pop.get_locations() < loc_2 )
        return 0;
    if ( pop.point_sub_pop(loc_1)->N() < 1 || pop.point_sub_pop(loc_2)->N() < 1 )
        return 0;
    int actual_sample_size = 0;
    if ( pop.point_sub_pop(loc_1)->N() < sample_size || pop.point_sub_pop(loc_2)->N() < sample_size )
        actual_sample_size = min (pop.point_sub_pop(loc_1)->N(), pop.point_sub_pop(loc_2)->N());
    else actual_sample_size = sample_size;

    pop.point_sub_pop(loc_1)->random_sample.clear();
    pop.point_sub_pop(loc_2)->random_sample.clear();
    pop.point_sub_pop(loc_1)->produce_random_sample(actual_sample_size);
    pop.point_sub_pop(loc_2)->produce_random_sample(actual_sample_size);

    string *str;
    str = new string [actual_sample_size];

    for (int i = 0; i < actual_sample_size/2; i ++)
    {
        str[i] = pop.point_sub_pop(loc_1)->get_genotype_string(pop.point_sub_pop(loc_1)->random_sample[i]);

    }
    for (int i = actual_sample_size/2; i < actual_sample_size; i ++)
    {
        str[i] = pop.point_sub_pop(loc_2)->get_genotype_string(pop.point_sub_pop(loc_1)->random_sample[i]);

    }
    int gen_size = str[0].size();

    for (int gen_pos = 0; gen_pos < gen_size; gen_pos ++)
    {
        int gen_no = 0;
        int k = 0;
        while (k == 0 && gen_no < actual_sample_size)
        {
            if (str[0][gen_pos] != str[gen_no][gen_pos])
            {
                result += 1;
                k = 1;
            }
            gen_no ++;
        }

    }
    return result;
}


struct genotype_stat
{


public:

    double average;
    double variance;


    genotype_stat(multi_population pop, int loc, int sample_size)
    {
        get_genotype_av(pop, loc, sample_size);

    }


private:

    double get_genotype_av(multi_population pop, int loc, int sample_size)
    {
        if (pop.get_locations() < loc || loc < 0)
        {
            average = variance = 0;
            return 0;
        }
        if (pop.point_sub_pop(loc)->N() == 0)
        {
            average = variance = 0;
            return 0;
        }


        string str;
        double aAverage = 0;
        double aVariance = 0;
        int normal = 1e-15;

        int gen_size = pop.point_sub_pop(loc)->L();
        int P[gen_size + 1];
        for (int i = 0; i < gen_size + 1; i ++)
        {
            P[i] = 0;
        }
        int actual_sample_size = 0;

        if ( pop.point_sub_pop(loc)->N() < 2 * sample_size )
            actual_sample_size = (pop.point_sub_pop(loc)->N() / 2);
        else actual_sample_size = sample_size;
        pop.point_sub_pop(loc)->random_sample.clear();
        pop.point_sub_pop(loc)->produce_random_sample(actual_sample_size);

        for (int i = 0; i < actual_sample_size; i ++)
        {
            int gen_sum = 0;
            str = pop.point_sub_pop(loc)->get_genotype_string(pop.point_sub_pop(loc)->random_sample[i]);

            for (int gen_pos = 0; gen_pos < gen_size; gen_pos ++)
            {
                if (str[gen_pos] == '1')
                    gen_sum ++;  // str[gen_pos];

            }

            P[gen_sum] ++;

        }

        for (int i=0; i < gen_size; i ++ )
        {
            cout << P[i] << " ";

        }
        cout << endl;


        for (int i = 0; i < gen_size + 1; i ++)
        {
            aAverage += i * P[i];
            aVariance += i * i * P[i];
            normal += P[i];
        }

        cout << actual_sample_size << endl;
        cout << normal << endl;

        average = aAverage/normal;
        variance = aVariance/normal - average * average;
        cout << average << "    " << variance << endl;
    }


};





/* MAIN */
int main(int argc, char  *argv[]){


        // ############################################
        //
        // Read the input parameters of the simulation
        //
        // ############################################
        double SELECTIVE = atof(argv[1]); //MIGRATION_RATE = atof(argv[1]);
        double POPULATION_SIZE = atof(argv[2]);
        const char * fileDir =  argv[3];
        int locations = atof(argv[4]);
       // int _L = atoi(argv[5]);
        // ############################################



        // ############################################
        //
        // Initialize the essential parameters for the simulation and load the population
        //
        // ############################################
        int Av_Num = 1;
        int generation = 0;
        int N = 20, N_pop = POPULATION_SIZE;
        int GENERATIONS = 1000;
        int SAMPLE_SIZE = N_pop / 10;
       /* if (MIGRATION_RATE != 0)
            GENERATIONS = 5 * max(double(POPULATION_SIZE), 1/MIGRATION_RATE);
        else
            GENERATIONS = POPULATION_SIZE;
*/
        double MIGRATION_RATE = 0.001;   //1;//SELECTIVE = 0.5/POPULATION_SIZE;



        int GENERATIONS_Eq = 10000; //100 * GENERATIONS;
        int n_o_traits = NUMBER_OF_TAITS;

        multi_population pop(locations, _L, n_o_traits);
        pop.set_migration_rate(MIGRATION_RATE);


        for (int i = 0; i < locations; i ++)
        {
            pop.point_sub_pop(i)->set_mutation_rate(MUTATION_RATE);
            pop.point_sub_pop(i)->outcrossing_rate = 0.0;
            pop.point_sub_pop(i)->crossover_rate = 0.0;//1e-2;
            pop.point_sub_pop(i)->recombination_model = FREE_RECOMBINATION;
        }
        vector <int> gen_loci;
        gen_loci.push_back(100);
        pop.track_locus_genealogy(gen_loci);
        // ############################################




        //#############################################
        //#
        //# Set the population and track the genealogy for the population
        //#
        //#############################################
        pop.point_sub_pop(0)->set_wildtype(N);		// start with a population of the right size
        for (int cur_loc=0; cur_loc < locations; cur_loc ++ ){
            pop.point_sub_pop(cur_loc)->carrying_capacity = N_pop;
        }
        pop.set_global_generation(generation);
        pop.submit_pop_genealogy();
        pop.genealogy.add_generation(pop.max_fitness());
        // ############################################





        //#############################################
        //#
        //# Set genotype -> trait -> fitness maps
        //#
        //#############################################
        vector <int> loci;

        for (int i = 0; i < locations; i ++)
        {


            //for (int trait_No = 0; trait_No < n_o_traits; trait_No ++)
            {
                for(int j=0; j < _L / 2; j++)
                {
                    loci.assign(1, j);
                    pop.point_sub_pop(i)->add_trait_coefficient(SELECTIVE / 2, loci, 0);
                    loci.clear();

                }

                for(int j=_L / 2; j < _L; j++)
                {
                    loci.assign(1, j);
                    pop.point_sub_pop(i)->add_trait_coefficient(0, loci, 0);
                    loci.clear();

                }


                for(int j=0; j < _L / 2; j++)
                {
                    loci.assign(1, j);
                    pop.point_sub_pop(i)->add_trait_coefficient(0, loci, 1);
                    loci.clear();

                }

                for(int j=_L / 2; j < _L; j++)
                {
                    loci.assign(1, j);
                    pop.point_sub_pop(i)->add_trait_coefficient(SELECTIVE / 2, loci, 1);
                    loci.clear();

                }



            }

            /*makes no sense in the temp calculations*/
            double * weights = new double[n_o_traits];
            for (int weight_No = 0; weight_No < n_o_traits; weight_No ++)
            {
                weights[weight_No] = 1;
            }
            pop.point_sub_pop(i)->set_trait_weights(weights);
            delete [] weights;
            /*=======================================*/

            pop.point_sub_pop(i)->update_traits();
            pop.point_sub_pop(i)->update_fitness();

        }

        stat_t fitstat;
        gsl_histogram *SFS = gsl_histogram_alloc(20);
        gsl_histogram_set_ranges_uniform(SFS,0,1);
        // ############################################




        // ##############################################
        // ##############################################
        // ##                                          ##
        // ##               EVOLVE CYCLE               ##
        // ##                                          ##
        // ##############################################
        // ##############################################
        for (int i = 0; i < GENERATIONS_Eq; i ++)
        {

            for (int cur_loc = 0 ; cur_loc < locations; cur_loc ++)
            {
                if (pop.point_sub_pop(cur_loc)->population_size > 0)
                {
                    pop.evolve(cur_loc, 1);
                }
            }

            pop.migrate();
            generation ++;
            pop.set_global_generation(generation);
            if (i % 10 == 0)
            {
                cout << "generation: " << generation << endl;
            }

            pop.submit_pop_genealogy();
            pop.genealogy.add_generation(pop.max_fitness());
        }
        // ##############################################
        // ##############################################



        //############################################################
        //#
        //# Take a sample from the population and output to the file
        //#
        //############################################################



        //############################################################
        //#
        //#  Create an output file for measured parameters
        //#  and set-up the format and comments for the output
        //#
        //############################################################
        stringstream  temp1 (stringstream::in | stringstream::out);
        stringstream  temp2 (stringstream::in | stringstream::out);
        temp1 << N_pop;
        temp2 << SELECTIVE;//MIGRATION_RATE;
        string tmp1 = temp1.str();
        string tmp2 = temp2.str();
        string outfilename;
        string outfileDir;
        outfileDir.append(fileDir).append("Meas");
        outfilename.append(outfileDir).append("_PopSize_").append(tmp1).append("_SelectionRate_").append(tmp2);

        const char * filename = outfilename.c_str();
        ofstream myfile;
        myfile.open (filename);
        myfile << "#Mutation rate:  " << MUTATION_RATE << endl << " #Population_size   " << POPULATION_SIZE <<  endl<<"#Sampling distance   " << GENERATIONS << endl;
        myfile << "#Number of migration events: "<< pop.number_of_migration_events << endl;
        myfile << "#Pairwise divergence within populations:" << endl;
        myfile << "#Location  1:   " << "Location  2: " << "Loc1 x Loc2:" << endl;
        myfile.close();
        //############################################################




        for (int SampleNum = 0; SampleNum < Av_Num; SampleNum++)
        {
            cout << "cycle " << SampleNum << endl;
            if (SampleNum != 0 )
            {


        // ##############################################
        // ##############################################
        // ##                                          ##
        // ##         EVOLVE CYCLE -- SAMPLING         ##
        // ##                                          ##
        // ##############################################
        // ##############################################


                for (int i = 0; i < GENERATIONS; i ++)
                {

                    for (int cur_loc = 0 ; cur_loc < locations; cur_loc ++)
                    {
                        if (pop.point_sub_pop(cur_loc)->population_size > 0)
                        {
                            pop.evolve(cur_loc, 1);
                        }
                    }

                    pop.migrate();
                    generation ++;

                    pop.set_global_generation(generation);
                    pop.submit_pop_genealogy();
                    pop.genealogy.add_generation(pop.max_fitness());
                }

            }
        // ##############################################
        // ##############################################


            ofstream outfile;
            const char * filename = outfilename.c_str();
            outfile.open(filename, std::ios_base::app);
            cout << filename << endl;
            //outfile << MIGRATION_RATE << " ";
         /*   for (int i = 0; i < locations; i ++)
            {
                outfile <<  pop.point_sub_pop(i)->get_pairwise_divergence(SAMPLE_SIZE) << "    ";
            }
            for (int i = 0; i < locations - 1; i ++)
            {
                for (int j = i + 1; j < locations; j ++)
                {
                    outfile <<  get_inter_pop_divergence (pop, i, j, SAMPLE_SIZE) << "  ";
                }
            }
            for (int i = 0; i < locations; i ++)
            {
                outfile <<  pop.point_sub_pop(i)->get_segregating_sites_num(SAMPLE_SIZE) << "    ";
            }

            for(int loc_1 = 0; loc_1 < locations; loc_1 ++)
            {
                for (int loc_2 = loc_1 + 1; loc_2 < locations; loc_2 ++)
                {
                    outfile << get_segregating_sites_num (pop, loc_1, loc_2, SAMPLE_SIZE) << "  ";
                }



            }*/

            outfile << endl;

            for (int locNum = 0; locNum < locations; locNum ++)
            {
                genotype_stat* gen_stat = new genotype_stat(pop, locNum, SAMPLE_SIZE);
                outfile << gen_stat->average << "   "<< gen_stat->variance << "    ";
                //outfile << gen_stat << endl;
                delete  gen_stat;
    //            outfile << pop.point_sub_pop(loc_1)->population_size << "   " << pop.point_sub_pop(loc_1)->get_fitness_statistics().mean << "   "<< pop.point_sub_pop(loc_1)->get_fitness_statistics().variance<< " ";
            }


    //        outfile << endl;
            outfile.close();

            stringstream  temp3 (stringstream::in | stringstream::out);
            stringstream  temp4 (stringstream::in | stringstream::out);
            temp3 << N_pop;
            temp4 << SELECTIVE;
            string tmp3 = temp3.str();
            string tmp4 = temp4.str();
            string out_tree_filename;
            string outfileDir_tree;
            outfileDir_tree.append(fileDir).append("Tree");
            out_tree_filename.append(outfileDir_tree).append("_PopSize_").append(tmp3).append("_SelectionRate_").append(tmp4);
            const char * tree_filename = out_tree_filename.c_str();

            ofstream myfile_tree;
            myfile_tree.open(tree_filename);
             for (unsigned int genlocus = 0; genlocus < gen_loci.size(); genlocus ++){
                //pop.genealogy.trees[genlocus].check_tree_integrity();
                myfile_tree <<"#PRINT ENTIRE TREE"<<endl;
                myfile_tree << pop.genealogy.trees[genlocus].print_newick()<<endl;

            }
             myfile_tree.close();



        }





return 0;
}

