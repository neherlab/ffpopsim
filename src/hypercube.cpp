/*
 *  hypercube_function.cpp
 *
 *  Created by Richard Neher on 01/27/09.
 */

#include "popgen.h"

//default constructor
hypercube::hypercube()
{
	mem=false;
	if (HC_VERBOSE) cerr<<"hypercube::hypercube(): constructing...!\n";
}

//constructor
hypercube::hypercube(int dim_in, int s)
{
	if (HC_VERBOSE) cerr<<"hypercube::hypercube(): constructing...!\n";
	set_up(dim_in);
}

//check consistency of input dimension and call allocation routine, seed for the rng is provided (s)
int hypercube::set_up(int dim_in, int s)
{
	if (dim_in<=8*sizeof(int) && dim_in>0)
	{
		if (HC_VERBOSE) cerr<<"hypercube::set_up(): setting up...!";
		dim=dim_in;
		mem=false;
		if (s==0) seed=time(NULL);
		else seed=s;
		if (HC_VERBOSE) cerr<<"done.\n";
		return allocate_mem();
	}
	else
	{
		cerr <<"hypercube: got: "<<dim_in<<", need positive and <32 dimension!\n";
		return HC_BADARG;
	}
}

//destructor
hypercube::~hypercube()
{
	if (HC_VERBOSE) cerr<<"hypercube::~hypercube(): destructing...!\n";
	if (mem) free_mem();
}

//allocate memory
int hypercube::allocate_mem()
{
	if (HC_VERBOSE) cerr<<"hypercube::allocate_mem(): allocating memory...";
	if (mem)
	{
		cerr <<"hypercube::allocate_mem(): memory already allocated, freeing and reallocating ...!\n";
		free_mem();
	}

	coeff=new double [1<<dim];  //(1<<dim)==2^dim
	order=new int [1<<dim];
	func=new double [1<<dim];

	if (func==NULL or coeff==NULL or order==NULL) //if allocation goes wrong, return error
	{
		cerr <<"hypercube::allocate_mem(): cannot allocate memory...!\n";
		free_mem();
	}

	//allocate and seed the random number generator, RNG is the rng-type
	rng=gsl_rng_alloc(RNG);
	gsl_rng_set(rng, seed);
	calc_order();
	mem=true;
	reset();
	if (HC_VERBOSE) cerr<<"done.\n";
	return 0;
}

//free memory
int hypercube::free_mem()
{
	if (!mem)
	{
		cerr <<"hypercube::free_mem(): no memory allocated...!\n";
		return 0;
	}
	delete [] coeff;
	delete [] func;
	delete [] order;
	mem=false;
	return 0;
}



/********** HC_FUNCTIONS THAT SET OR INCREMENT THE HC_COEFFICIENTS *******/
//initialize the coefficients such that terms of order k contribute a variance var[k]
//coefficients are gaussian random numbers.
int hypercube::gaussian_coefficients(double* var, bool add)
{
	int temp;
	double sigma;

	if (add and state==HC_FUNC) fft_func_to_coeff();

	coeff[0]=0;		//constant coefficient is set to zero
	if (HC_VERBOSE) cerr<<"hypercube::init_gauss_var(): initialize with specified variance...";

	//loop over all coefficients
	for (int k=1; k<(1<<dim); k++)
	{
		sigma=sqrt(var[k]);
		temp=order[k];
		if (var[temp]>0)
		{
			//init with gaussian random number of zero mean and unit  variance
			if (add) coeff[k]+=(gsl_ran_gaussian(rng, sqrt(var[temp])));
			else coeff[k]=(gsl_ran_gaussian(rng, sqrt(var[temp])));
		}else
		{
			if (!add) coeff[k]=0.0;
		}
	}
	fft_coeff_to_func();
	if (HC_VERBOSE) cerr<<"done.\n";
	return 0;
}

int hypercube::additive(double* additive_effects, bool add)
{
	if (add and state==HC_FUNC) fft_func_to_coeff();

	if (HC_VERBOSE) cerr<<"hypercube::additive(): initialize with specified variance...";
	if (!add) reset();
	//loop over all coefficients
	for (int locus=0; locus<dim; locus++)
	{
		if (add) coeff[1<<locus]+=additive_effects[locus];
		else coeff[1<<locus]=additive_effects[locus];
	}
	fft_coeff_to_func();
	if (HC_VERBOSE) cerr<<"done.\n";
	return 0;
}


/*
 * Functions that return the maximum of the function or the corresponding genotype
 */
int hypercube::argmax()
{
	if (state==HC_COEFF) fft_coeff_to_func();

	int max_index=0;
	double max_value=func[0];
	for (int i=0; i<(1<<dim); i++){
		if (func[i]>max_value) {max_index=i; max_value=func[i];}
	}
	return max_index;
}

double hypercube::valuemax()
{
	if (state==HC_COEFF) fft_coeff_to_func();
	int max_index=0;
	double max_value=func[0];
	for (int i=0; i<(1<<dim); i++){
		if (func[i]>max_value) {max_index=i; max_value=func[i];}
	}
	return max_value;
}


/******* HC_FUNCTIONS THAT INITIALIZE OR INCREMENT THE HC_FUNCTION ON 2^L*******/

//init with constant spectrum, gaussian random numbers
int hypercube::init_rand_gauss(double sigma, bool add)
{
	if (add and state==HC_COEFF) fft_coeff_to_func();
	if (add)  {for (int i=0; i<(1<<dim); i++){ func[i]+=gsl_ran_gaussian(rng,sigma);}}
	else {for (int i=0; i<(1<<dim); i++){ func[i]=gsl_ran_gaussian(rng,sigma);}}

	return fft_func_to_coeff();
}

int hypercube::init_list(vector <index_value_pair> iv, bool add){
	if (add==false){
		reset();
	}else if (state==HC_COEFF)	fft_coeff_to_func();

	for (unsigned int pair=0; pair<iv.size(); pair++){
		func[iv[pair].index]=iv[pair].val;
	}
	return fft_func_to_coeff();
}


int hypercube::init_coeff_list(vector <index_value_pair> iv, bool add){
	if (add==false){
		reset();
	}else if (state==HC_FUNC)	fft_func_to_coeff();

	for (unsigned int pair=0; pair<iv.size(); pair++){
		coeff[iv[pair].index]=iv[pair].val;
	}
	return fft_coeff_to_func();
}


void hypercube::calc_order()
{
	int spin;

	order[0]=0;
	if (HC_VERBOSE) cerr<<"hypercube::calc_order(): make a list of the order of different coefficients...";

	//loop over all coefficients
	spin=-1;	//auxilliary variable
	for (int k=1; k<(1<<dim); k++)
	{
		if (k==(1<<(spin+1))) spin++;
		order[k]=1+order[k-(1<<spin)];	//the order of coefficient k is 1+(the order of coefficient[k-2^spin])
	}
}
/*
 * UTILITY HC_FUNCTIONS
 */
int hypercube::normalize(double targetnorm){
	double norm=0;
	for (int i=0; i<(1<<dim); i++){
		norm+=func[i];
	}
	return scale(targetnorm/norm);
}

//rescale the coefficients
int hypercube::scale(double scale)
{
	if (!mem)
	{
		cerr <<"hypercube::scale: allocate memory first!\n";
	}
	if (state==HC_FUNC) fft_func_to_coeff();
	else if (state==HC_COEFF) fft_coeff_to_func();

	for(int i=0; i<(1<<dim); i++)
	{
		coeff[i]*=scale;
		func[i]*=scale;
	}
	return 0;
}
//rescale the function
int hypercube::shift(double shift)
{
	if (!mem)
	{
		cerr <<"hypercube::shift: allocate memory first!\n";
	}
	if (state==HC_FUNC) fft_func_to_coeff();
	else if (state==HC_COEFF) fft_coeff_to_func();

	coeff[0]+=shift;
	for(int i=0; i<(1<<dim); i++)
	{
		func[i]+=shift;
	}
	return 0;
}

//reset function values and coefficients
int hypercube::reset()
{
	if (!mem)
	{
		cerr <<"hypercube::reset: allocate memory first!\n";
	}
	for(int i=0; i<(1<<dim); i++)
	{
		coeff[i]=0;
		func[i]=0;
	}
	state=HC_FUNC_EQ_COEFF;
	return 0;
}


/******** FFT and its INVERSE **********/

//perform the fourier transform to calculate the function from the coefficients
int hypercube::fft_coeff_to_func()
{
	if (HC_VERBOSE) cerr<<"hypercube::fft_coeff_to_func(): calculate function from coefficients....";
	int i,k, spinmask;
	double *temp=new double [(1<<dim)];
	double *swap_pointer;

	if (!mem)
	{
		cerr <<"hypercube::fft_coeff_to_func(): allocate memory first!\n";
		return HC_MEMERR;
	}

	//copy coefficients into the memory hold for the function values
	for (i=0; i<(1<<dim); i++)
	{
		func[i]=coeff[i];
	}

	//loop over all spins
	for (k=0; k<dim; k++)
	{
		spinmask=1<<k;	// integer value with a 1 at position k in binary representation
		for (i=0; i<(1<<dim); i++)	//loop over the hypercube
		{
			if (i&spinmask)	//if spin k is set
			{
				//pull out what differs between (spin k)=0 (at i-spinmask) and (spin k)=1
				temp[i]=(func[i-spinmask]+func[i]);
			}
			else			//if not
			{
				//pull out what is independent of spin k
				temp[i]=(func[i]-func[i+spinmask]);
			}
		}
		//swap func and temp
		swap_pointer=func;
		func=temp;
		temp=swap_pointer;
	}
	delete [] temp;
	if (HC_VERBOSE) cerr<<"done!\n";

	state=HC_FUNC_EQ_COEFF;	//set state to
	return 0;
}

//perform the inverse fourier transform to calcuate the coefficients from the function
int hypercube::fft_func_to_coeff()
{
	int i,k, spinmask;
	double *temp=new double [(1<<dim)];
	double *swap_pointer;
	if (!mem)
	{
		cerr <<"hypercube::fft_func_to_coeff(): allocate memory first!\n";
		return HC_MEMERR;
	}


	//copy coefficients into the memory hold for the function values
	for (i=0; i<(1<<dim); i++)
	{
		coeff[i]=func[i];
	}

	//loop over all spins
	for (k=0; k<dim; k++)
	{
		spinmask=1<<k;	// integer value with a 1 at position k in binary representation
		for (i=0; i<(1<<dim); i++)	//loop over the hypercube
		{
			if (i&spinmask)	//if spin k is set
			{
				//pull out what differs between (spin k)=0 (at i-spinmask) and (spin k)=1
				temp[i]=0.5*(coeff[i]-coeff[i-spinmask]);
			}
			else			//if not
			{
				//pull out what is independent of spin k
				temp[i]=0.5*(coeff[i+spinmask]+coeff[i]);
			}
		}
		//swap func and temp
		swap_pointer=coeff;
		coeff=temp;
		temp=swap_pointer;
	}
	delete [] temp;
	return 0;
}


/******** END FFT ************/

/********* INPUT OUTPUT OF HC_FUNCTION AND HC_COEFFICIENTS **********/
//read coefficients from a stream
int hypercube::read_coeff(istream &in)
{
	int i;
	if (in.bad())
	{
		cerr <<"hypercube::read_coeff: bad stream\n";
		return HC_BADARG;
	}
	i=0;
	while(in.good() && i<(1<<dim))
	{
		in >>coeff[i];
		i++;
	}
	if (i<(1<<dim))
	{
		cerr <<"hypercube::read_coeff: file end reached after "<<i<<" values!\n";
		return HC_BADARG;
	}
	return 0;

}

//write func to a stream
int hypercube::write_func(ostream &out)
{
	int i;
	if (out.bad())
	{
		cerr <<"hypercube::write_func: bad stream\n";
		return HC_BADARG;
	}
	i=0;
	while(out.good() && i<(1<<dim))
	{
		out <<func[i]<<endl;
		i++;
	}
	if (i<(1<<dim))
	{
		cerr <<"hypercube::write_func: error while writing!\n";
		return HC_BADARG;
	}
	return 0;

}

//read coefficients from a stream
int hypercube::read_func(istream &in)
{
	int i;
	if (in.bad())
	{
		cerr <<"hypercube::read_func: bad stream\n";
		return HC_BADARG;
	}
	i=0;
	while(in.good() && i<(1<<dim))
	{
		in >>func[i];
		i++;
	}
	if (i<(1<<dim))
	{
		cerr <<"hypercube::read_func: file end reached after "<<i<<" values!\n";
		return HC_BADARG;
	}
	return 0;
}

//read coefficients from a stream
int hypercube::read_func_labeled(istream &in)
{
	int i, count;
	char gt[dim+2];
	if (in.bad())
	{
		cerr <<"hypercube::read_func_labeled: bad stream\n";
		return HC_BADARG;
	}
	count=0;
	while(in.good() && count<(1<<dim))
	{
		i=0;
		in >>gt;
		for (int k=0; k<dim; k++)
			if (gt[k]=='1') i+=1<<k;
		in >>func[i];
		count++;
		in.ignore(100, '\n');
	}
	if (count<(1<<dim))
	{
		cerr <<"hypercube::read_func: file end reached after "<<i<<" values!\n";
		return HC_BADARG;
	}
	return 0;
}


//write coeff to a stream
int hypercube::write_coeff(ostream &out, bool label)
{
	int i,k;
	if (out.bad())
	{
		cerr <<"hypercube::write_coeff: bad stream\n";
		return HC_BADARG;
	}
	i=0;
	while(out.good() && i<(1<<dim))
	{
		if (label)
		{
			for (k=0; k<dim; k++)
			{
				if (i&(1<<k)) out <<'1';
				else out <<'0';
			}
			out <<"  ";
		}
		out <<coeff[i]<<endl;
		i++;
	}
	if (i<(1<<dim))
	{
		cerr <<"hypercube::write_coeff: error while writing!\n";
		return HC_BADARG;
	}
	return 0;

}

/********* END OF INPUT OUTPUT OF HC_FUNCTION AND HC_COEFFICIENTS **********/


//function to test the fourer transform, output is written to the error stream
int hypercube::test()
{
	int tp, sign;
	double res;
	cerr <<"hypercube::test()...\n";
	for (int i=0; i<100; i++)
	{
		res=0;
		//draw a random integer
		tp=gsl_rng_uniform_int(rng,1<<dim);
		for (int c=0; c<(1<<dim); c++)		//assemble the function value from its spins, go over all coeffcients, index c
		{
			sign=1;
			for (int k=0; k<dim; k++)
			{
				if (!(tp&(1<<k)) && (c&(1<<k)))		//determine the sign the coeff enters
				{
					sign*=-1;						//multiply by minus one, every time a bit is not set although it the contribution depends on the bit
				}
			}
			res+=sign*coeff[c];						//add the coeff with the correct sign
		}
		cerr <<setw(15)<<func[tp]<<setw(15)<<res<<setw(15)<<func[tp]-res<<endl;	//output result, function value and difference, the latter should be zero
	}
	cerr <<"done\n";
	return 0;
}

