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

class environmental_hh : public haploid_highd{

    double sigma;
    virtual double trait_function(double aTrait1, double aTrait2){
        double t1 = aTrait1;
        double t2 = aTrait2;
        double rho = 0;
        rho = sqrt(pow((t1 + 1e-5), 2) + pow((t2 + 1e-5), 2));
        double d_phi = acos(((t1 + 1e-5) * cos(phi_0) + (t2 + 1e-8) * sin(phi_0)) / rho);
        return trait_weights[0]  * exp(1) *  rho / (rho_0 + 1e-8) * exp(- rho / (rho_0 + 1e-8)) * exp(- pow((d_phi), 2) / sigma);
    }
    double rho0;
    double phi0;
public:

    //get/set methods for fitness map parameters
    void setPhi0 (double new_phi0){
        double new_phi = new_phi0;
        while (new_phi < 0 ){
                new_phi = new_phi + 2 * 3.14;
        }
        while (new_phi > 2 * 3.14){
                new_phi = new_phi - 2 * 3.14;
        }
        phi_0 = new_phi;
        return;
    }
    double get_phi_0 (){return phi_0;}
    void setRho0(double new_rho0){rho0 = new_rho0; return;}
    double getRho0(){return rho0;}
    //This method doeas the job of trait -> fitness conversion
    virtual void calc_individual_fitness_from_traits(clone_t &tempgt) {tempgt.fitness = trait_function(tempgt.trait[0], tempgt.trait[1]);}
    environmental_hh(){
        phi0 = 3.14 / 4;
        rho0 = 8;
        sigma = 0.1;
    }
    ~environmental_hh(){}
};


class multi_population
{
    double MAX_MIGRATION_RATE;
    double critical_migration_rate;
    //vector <haploid_highd> sub_population;
    vector <environmental_hh*> sub_population;
    //vector<int> population_sizes;
    //vector <haploid_highd> sub_population;
    int number_of_loci;
    int number_of_traits;
    int number_of_locations;
    int track_genealogy;
    double fitness_max;

    //double migration_rate;
    vector < vector < double > > migration_rates;
    vector < double > mutation_rates;
    vector< int > carrying_capacities;

    void add_migrating_clone_to_genealogy(int locusIndex, int old_location, int new_location,  int dest, int parent, int left, int right, int cs, int n);
    void set_migration_rates(vector<double> new_migration_rates, int source);
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
    void update_population_size();

public:

    void reset();
    environmental_hh * point_sub_pop(int i){return sub_population[i];}
    multi_population(int new_locations, int L_in, int n_o_traits = 1, int rng_seed = 0);
    ~multi_population();

    //parameteres
    int get_number_of_locations () {return number_of_locations;}
    int N(int i=-1);
    int L(){return number_of_loci;}
    int get_generation(){return generation;}
    double max_fitness();
    int number_of_migration_events;

    //migration
    int set_migration_rates(vector< vector<double> > new_migration_rates);
    vector< vector<double> > get_migration_rates() {return migration_rates;}
    double get_critical_migration_rate(){return critical_migration_rate;}
    void set_critical_migration_rate(double new_rate);


    //mutation
    void set_mutation_rates(vector<double> rates);
    vector<double> get_mutation_rates(){return mutation_rates;}
    double get_mutation_rate(int location);

    //carrying capacities
    void set_carrying_capacities(vector<int> capacities);
    double get_carrying_capacity(int location);
    vector<int> get_carrying_capacities(){return carrying_capacities;}
    void set_tree_sample(vector <int> tree_sample){for (size_t i=0; i<number_of_locations; i++) sub_population[i]->tree_sample=tree_sample[i];}

    //trais and phenotype
    int get_number_of_traits(){return number_of_traits;}
    void set_environment(double rho, double phi, int location=-1);
    void set_trait_coefficient(double coefficient, vector<int> loci, int trait_no);
    void set_trait_weights(double* weights);
    double get_trait_weight(int trait_no){return point_sub_pop(0)->get_trait_weight(trait_no);}
    void update_fitness();

    void set_outcrossing_rate(double o_rate);
    double get_outcrossing_rate(){return point_sub_pop(0)->outcrossing_rate;}

    void set_crossover_rate(double c_rate);
    double get_crossover_rate(){return point_sub_pop(0)->crossover_rate;}

    void set_recombination_model(int r_model);
    int  get_recombination_model(){return point_sub_pop(0)->recombination_model;}

    //genealogy
    multi_locus_genealogy genealogy;
    void track_locus_genealogy(vector<int> loci);

    //get data
    boost::dynamic_bitset<> get_genotype(clone_t clone);
    clone_t get_random_clone(int location = -1);

    //evolution
    int evolve(int gen = 1);

    //init population
    void set_random_genotype(int N_in);
    void set_wildtype(int N_in);
protected:

    int determine_number_of_migrants(haploid_highd sub_population);
    int pickup_migrant(haploid_highd sub_population);
    int transfer_clone(int sub_pop_source, int sub_pop_destination, int source);
    void migrate();
    void migrate(int source);
};


#endif // MULTI_POPULATION_H
