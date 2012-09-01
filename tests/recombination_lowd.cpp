/**
 * @file recombination_lowd.cpp
 * @brief Test recombination routine for haploid_lowd
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-08-28
 *
 * Tests for the recombination routine in haploid_lowd are performed,
 * for the general case and for single crossover.
 */
/* Include directives */
#include <string>
#include "ffpopsim_lowd.h"
#define LOWD_BADARG -1354341
#define NOTHING 1e-10

/* Be verbose? */
#define LOWD_VERBOSE 1

/* Declaration */
class haploid_lowd_test : public haploid_lowd {
public:
	// construction / destruction
	haploid_lowd_test(int L=1, int rngseed=0) : haploid_lowd(L, rngseed){};
	~haploid_lowd_test();

	//testing
	int test_recombinant_distribution();
	int test_recombination(double *rec_rates);
	int mutation_drift_equilibrium(double** mutrates);
	int test_single_crossover_set_rates();
};


/* Implementation */
/**
 * @brief Default destructor
 */
haploid_lowd_test::~haploid_lowd_test() {}

/**
 * @brief Test the recombination routine using Fourier transforms
 *
 * @returns zero if both routines agree, -1 otherwise
 *
 * Debugging routine: calculates the distribution of recombinants explicitly and
 * compares the result to the recombinant distribution obtained via fourier transform
 */
int haploid_lowd_test::test_recombinant_distribution(){
	double *test_rec;
	double dev=0;
	//allocate memory for the recombinant distribution calculated step-by-step
	test_rec=new double [(1<<number_of_loci)];
	int mother, father;
	//now calculate the recombinant distribution from pairs of parents.
	int gt1, gt2, rec_pattern, inherited;

	cout <<"test_recombinant_distribution(): start..."<<endl;

	if (recombination_model == FREE_RECOMBINATION){
		//calculate recombinants the efficient way
		calculate_recombinants_free();
		for (gt1=0; gt1<(1<<number_of_loci); gt1++){	//target genotype
			test_rec[gt1]=0.0;							//initialize
			//loop over all recombination patterns (equal probability)
			for (rec_pattern=0; rec_pattern<(1<<number_of_loci); rec_pattern++){
				//loop over the parts of the maternal and paternal genomes not inherited
				for (gt2=0; gt2<(1<<number_of_loci); gt2++){
					//construct maternal and paternal genotypes
					mother=(gt1&(rec_pattern))+(gt2&(~rec_pattern));
					father=(gt1&(~rec_pattern))+(gt2&(rec_pattern));
					//increment the rec distribution
					test_rec[gt1]+=population.func[mother]*population.func[father];
				}
			}
			//normalize
			test_rec[gt1]*=1.0/(1<<number_of_loci);
			cout <<gt1<<"  "<<test_rec[gt1]<<"  "<<recombinants.func[gt1]<<endl;
			//sum up all deviations
			dev+=(test_rec[gt1]-recombinants.func[gt1])*(test_rec[gt1]-recombinants.func[gt1]);
		}
	} else if (recombination_model == SINGLE_CROSSOVER) {
		//calculate recombinants the efficient way
		calculate_recombinants_single();
		for (gt1=0; gt1<(1<<number_of_loci); gt1++) test_rec[gt1]=0.0;

		// the recombination patterns are L-1, plus one "nothing happens" pattern
		for (int crossover_point=0; crossover_point < number_of_loci; crossover_point++){

			// set the recombination pattern
			rec_pattern = 0;
			for (int locus=0; locus != number_of_loci - 1 - crossover_point; locus++)
				rec_pattern += (1<<(number_of_loci - 1 - locus));

			// include both patterns, 00111 and 11000
			for (int first_digit=0; first_digit <= 1; first_digit++) {
				if(first_digit)	rec_pattern = (1<<number_of_loci) - 1 - rec_pattern;

				// print debugging info
				for (int locus=number_of_loci - 1; locus >= 0; locus--) {
					cout<<((rec_pattern&(1<<locus))?1:0);	
				}
				cout<<"\t"<<rec_pattern<<endl;
				for (int locus=0; locus<number_of_loci-1; locus++) {
					cout<<(locus<=crossover_point?" ":(locus==(crossover_point+1)?"|":" "));					
				}
				cout<<endl<<endl;

				// iteratre over mothers and fathers
				for (mother=0; mother<(1<<number_of_loci); mother++){
					for (father=0; father<(1<<number_of_loci); father++){

						//contribution is weighted by the probability of this particular recombination pattern
						//this got calculated and stored in recombination_patterns[(1<<number_of_loci)-1]
						inherited = (mother&(rec_pattern)) + (father&(~rec_pattern));
//						if(inherited == 1)
//							cout<<crossover_point<<" "<<mother<<" "<<father<<endl;

						test_rec[inherited]+=recombination_patterns[0][crossover_point]*population.func[mother]*population.func[father];
					}
				}
			}
		}
		// normalize
		double sum=0;
		for (gt1=0; gt1<(1<<number_of_loci); gt1++) sum+=test_rec[gt1];
		for (gt1=0; gt1<(1<<number_of_loci); gt1++) test_rec[gt1]/=sum;
		// output
		cout<<"Genotype slow fast"<<endl;
		for (gt1=0; gt1<(1<<number_of_loci); gt1++) {
			cout <<gt1<<"  "<<test_rec[gt1]<<"  "<<recombinants.func[gt1]<<endl;
			dev+=(test_rec[gt1]-recombinants.func[gt1])*(test_rec[gt1]-recombinants.func[gt1]);
		}
	} else { 	//same as above, only the individual contribution
		//calculate recombinants the efficient way
		calculate_recombinants_general();
		for (gt1=0; gt1<(1<<number_of_loci); gt1++){
			test_rec[gt1]=0.0;
			for (rec_pattern=0; rec_pattern<(1<<number_of_loci); rec_pattern++){
				for (gt2=0; gt2<(1<<number_of_loci); gt2++){
					mother=(gt1&(rec_pattern))+(gt2&(~rec_pattern));
					father=(gt1&(~rec_pattern))+(gt2&(rec_pattern));
					//contribution is weighted by the probability of this particular recombination pattern
					//this got calculated and stored in recombination_patterns[(1<<number_of_loci)-1]
					test_rec[gt1]+=recombination_patterns[(1<<number_of_loci)-1][rec_pattern]*population.func[mother]*population.func[father];
				}
			}
			cout <<gt1<<"  "<<test_rec[gt1]<<"  "<<recombinants.func[gt1]<<endl;
			dev+=(test_rec[gt1]-recombinants.func[gt1])*(test_rec[gt1]-recombinants.func[gt1]);
		}
	}
	delete [] test_rec;
	if (dev>1e-9){
		cout <<"Deviation between explicit and fourier transform version! "<<dev<<endl;
		return -1;
	}else{
		cout <<"Explicit and fourier transform version agree to "<<dev<<endl;
		return 0;
	}
	return 0;
}

/**
 * @brief Test the recombination routine extensively
 *
 * @param rec_rates recombination rates used for testing
 *
 * @returns zero (but look at the stdout)
 *
 * Debugging routine: produces random genotypes configurations and test whether they recombine correctly.
 */
int haploid_lowd_test::test_recombination(double *rec_rates){

	//calculate the genetic map, i.e. cumulative recombination rates
	double* cumulative_rates=new double [number_of_loci+1];
	cumulative_rates[0]=0.0;
	for (int locus =1; locus<number_of_loci+1; locus++) cumulative_rates[locus]=cumulative_rates[locus-1]+rec_rates[locus-1];

	//initialize the internal recombination rates
	set_recombination_rates(rec_rates);

	//initialize the population randomly and test the recombination procedure
	population.set_state(HC_FUNC);
	for (int r=0; r<1; r++){
		for (int i=0; i<(1<<number_of_loci); i++){
			population.func[i]=gsl_rng_uniform(rng);
		}
		population.normalize();
		test_recombinant_distribution();
	}
	population.set_state(HC_FUNC);
	for (int i=0; i<(1<<number_of_loci); i++){
		population.func[i]=gsl_rng_uniform(rng);
	}
	population.normalize();

	//study the decay of cumulants from the randomly initialized initialized population
	//output header
	cout <<"\n\nRatio of the cumulants and the expected decay curve, should be constant. Last column shows dynamic range\n";
	cout <<"Generation  ";
	for (int l1=0; l1<number_of_loci; l1++){
		for(int l2=0; l2<l1; l2++){
			cout <<setw(13)<<l1<<" "<<l2;
		}
	}
	cout <<setw(15)<<"exp(-rmax*t)";
	cout<<'\n';
	//for a thousand time steps, recombine and watch the cumulants decay
	for (int g=0; g<1000; g++){
		if (g%100==0){ //output every hundred generations
			cout <<setw(10)<<g;
			for (int l1=0; l1<number_of_loci; l1++){
				for(int l2=0; l2<l1; l2++){
					cout <<setw(15)<<4 * get_LD(l1,l2)*exp(g*0.5*(1.0-exp(-2.0*(cumulative_rates[l1+1]-cumulative_rates[l2+1]))));
				}
			}
			cout <<setw(15)<<exp(-g*0.5*(1.0-exp(-2.0*(cumulative_rates[number_of_loci]-cumulative_rates[1]))));
			cout<<'\n';
		}
		recombine();
	}
	delete cumulative_rates;
	return 0;
}

/**
 * @brief Test the mutation-drift equilibrium with diffusion theory
 *
 * @param mu mutation rates
 *
 * @returns zero (but look at the stdout)
 */
int haploid_lowd_test::mutation_drift_equilibrium(double **mu){
	set_mutation_rates(mu);

	//init population and recombination rates
	double *af=new double[number_of_loci];;
	double *recrates=new double[number_of_loci];;
	for (int i=0; i<number_of_loci; i++){
		af[i]=0;
		recrates[i]=10;
	}
	set_allele_frequencies(af, 1000);
	//allocate histograms to store allele frequency distributions
	gsl_histogram **mutfreq=new gsl_histogram* [number_of_loci];
	for (int locus=0; locus<number_of_loci; locus++){
		mutfreq[locus]=gsl_histogram_alloc(100);
		gsl_histogram_set_ranges_uniform(mutfreq[locus], -1,1);
	}

	//equilibrate for 2N generations
	for (int gen=0; gen<2*carrying_capacity; gen++){
		mutate();
		resample();
	}
	//take 100000 samples every 1000 generations (assumes population is of order 1000)
	for (int r=0; r<100000; r++){
		for (int gen=0; gen<1000; gen++){
			mutate();
			resample();
		}

		for (int locus=0; locus<number_of_loci; locus++){
			gsl_histogram_increment(mutfreq[locus], get_chi(locus));
		}
	}

	//output: normalized histograms as well as theoretical expectation from diffusion theory
	//calculate norm of distributions first, output below.
	double upper, lower;
	double* histogramnorm=new double [number_of_loci];
	double* theorynorm=new double [number_of_loci];
	for (int locus=0; locus<number_of_loci; locus++){
		histogramnorm[locus]=0;
		theorynorm[locus]=0;
		for (int i=0; i<100; i++){
			gsl_histogram_get_range(mutfreq[locus], i, &lower, &upper);
			histogramnorm[locus]+=gsl_histogram_get(mutfreq[locus], i);
			theorynorm[locus]+=pow(0.5*(1+0.5*(upper+lower)), 2*carrying_capacity*mu[0][locus]-1)*pow(0.5*(1-0.5*(upper+lower)), 2*carrying_capacity*mu[1][locus]-1);
		}
	}
	for (int i=0; i<100; i++){
		gsl_histogram_get_range(mutfreq[0], i, &lower, &upper);
		cout <<setw(15)<<0.5*(upper+lower);
		for (int locus=0; locus<number_of_loci; locus++){
			cout <<setw(15)<<gsl_histogram_get(mutfreq[locus], i)/histogramnorm[locus]
					<<setw(15)<<pow(0.5*(1+0.5*(upper+lower)), 2*carrying_capacity*mu[0][locus]-1)*pow(0.5*(1-0.5*(upper+lower)), 2*carrying_capacity*mu[1][locus]-1)/theorynorm[locus];
		}
		cout <<endl;
	}
	return 0;
}

int haploid_lowd_test::test_single_crossover_set_rates() {
	cout<<"test_single_crossover_set_rates(): start..."<<endl;

	// print number of loci
	cout<<"L = "<<number_of_loci<<endl<<endl;

	// prepare rates
	cout<<"rates:";
	double *rec_rates = new double[number_of_loci - 1];
	for(int i=0; i < number_of_loci - 1; i++) {
		rec_rates[i] = 1e-2 / pow(10, (double)i);
		cout<<" "<<rec_rates[i];
	}
	cout<<endl<<endl;

	// set rates
	set_recombination_rates(rec_rates, SINGLE_CROSSOVER);
	//set_recombination_rates(rec_rates, CROSSOVERS);

	// check patterns
	vector <int> ii;
	int set_size;
	cout<<"Patterns\trates"<<endl;
	for (int subset=(1<<number_of_loci) - 1; subset != 0; subset--) {
		set_size = fitness.order[subset];
		ii.clear();
		for (int locus=0; locus < number_of_loci; locus++) {
			if ((subset&(1<<locus)))
				ii.push_back(locus);
		}

		for (int locus=0; locus < set_size; locus++) {
			for(int i=0; i < set_size; i++) cout<<((i <= locus)?0:1);
			cout<<" ";
			for(int i=0; i < set_size; i++) cout<<((i <= locus)?1:0);
			cout<<"\t"<<recombination_patterns[subset][locus]<<endl;
			for(int i=0; i < set_size; i++) cout<<ii[i];
				
			cout<<endl<<endl;
		}
	}
	cout<<"(NULL)\t"<<recombination_patterns[0][0]<<endl;

	cout <<"test_recombinant_distribution(): end."<<endl;

	// free memory
	delete [] rec_rates;
	return 0;

}


/* MAIN */
int main(int argc, char **argv){
	int status=0;
	if (argc > 1) {
		cout<<"Usage: "<<argv[0]<<endl;
		status = 1;
	} else {

		int L = 5;
		int N = 10000;
		haploid_lowd_test pop(L, 1);

		// start with wildtype and a single mutant
		index_value_pair_t ivp(0, N/4);
		vector<index_value_pair_t> gts;
		gts.push_back(ivp);
		ivp.index = 1;
		gts.push_back(ivp);
		ivp.index = 6;
		gts.push_back(ivp);
		ivp.index = 15;
		gts.push_back(ivp);
		pop.set_genotypes(gts);

		// test recombination routine
		//if (pop.test_single_crossover_set_rates()) status += 1;
		if (pop.test_recombinant_distribution()) status += 1;

	}
	cout<<"Number of errors: "<<status<<endl;
	return status;
}
