#include "ffpopsim_lowd.h"
int main() {
	haploid_lowd pop(4);
	pop.set_wildtype(1e10);
	pop.set_mutation_rates(1e-2);
	pop.evolve(1000);
	for (int i=0; i<(1<<4); i++)
		cout<<pop.get_genotype_frequency(i)<<endl;
	return 0;
}
