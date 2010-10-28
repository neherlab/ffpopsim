/*
 * sample.cpp
 *
 *  Created on: Aug 28, 2008
 *      Author: neher
 */

#include "popgen.h"

sample::sample() {
	mem_dis=false;
	mem_values=false;
	with_range=false;
	bins=0;
	number_of_values=0;
	values=NULL;
}

sample::~sample() {
	if (mem_dis)
		gsl_histogram_free(distribution);
	if (mem_values)
		delete [] values;
}

int sample::set_up(int n)
{
	if (n<1)
	{
		cerr <<"sample::set_up(): number of values has to be greater than zero! Got: "<<bins<<endl;
		return SAMPLE_ERROR;
	}
	else
	{
		number_of_values=n;
		values=new double [number_of_values];
		mem_values=true;
	}
	return 0;
}

int sample::set_distribution(int bins_in)
{
	if (bins_in<1)
	{
		cerr <<"sample::set_distribution(): number of bins has to be greater than zero! Got: "<<bins<<endl;
		return SAMPLE_ERROR;
	}
	else
	{
		bins=bins_in;
		distribution=gsl_histogram_alloc(bins);
	}
	return 0;
}

int sample::calc_mean()
{
	double v;
	if (mem_values==false)
	{
		cerr <<"sample::calc_mean(): Set values first!"<<endl;
		return SAMPLE_ERROR;
	}
	mean=0;
	for (int i = 0; i < number_of_values; ++i) {
		v=values[i];
		if (with_range)
		{
			v=(v>range_min)?v:range_min;
			v=(v<range_max)?v:range_max;
		}
		mean += v;
	}
	mean/=number_of_values;
	return 0;
}

int sample::calc_variance()
{
	double v;
	if (mem_values==false)
	{
		cerr <<"sample::calc_variance(): Set values first!"<<endl;
		return SAMPLE_ERROR;
	}
	calc_mean();
	variance=0;
	for (int i = 0; i < number_of_values; ++i) {
		v=values[i];
		if (with_range)
		{
			v=(v>range_min)?v:range_min;
			v=(v<range_max)?v:range_max;
		}
		variance += v*v;
	}
	variance/=(number_of_values-1);
	variance-=mean*mean*number_of_values/(number_of_values-1);
	return 0;
}

int sample::calc_distribution()
{
	if (mem_values==false)
	{
		cerr <<"sample::calc_distribution(): Set values first!"<<endl;
		return SAMPLE_ERROR;
	}
	double min=values[0], max=values[0];
	if (with_range)
	{
		min=range_min;
		max=range_max;
	}
	else
	{
	for (int i = 1; i < number_of_values; ++i) {
		if (values[i]<min) min=values[i];
		else if (values[i]>max) max=values[i];
	}
	}
	if (min==max) {max++; min--;}
	//set ranges and reset histogram bins to zero
	gsl_histogram_set_ranges_uniform(distribution, min, max);
	for (int i = 0; i < number_of_values; ++i) {
		gsl_histogram_increment(distribution, values[i]);
	}
	return 0;
}

//print allele frequencies to stream
int sample::print_distribution(ostream &out)
{
	if (out.bad())
	{
		cerr <<"sample::print_distribution(): bad stream\n";
		return SAMPLE_ERROR;
	}
	calc_distribution();
	for (int l=0; l<bins; l++) out <<setw(15)<<gsl_histogram_get(distribution, l);
	out <<endl;
	return 0;
}


