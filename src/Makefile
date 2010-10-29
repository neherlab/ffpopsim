CFLAGS= -I ~/include -O2 -c -g


PopGenLib:	haploid_gt_dis.o    haploid_population.o hypercube_function.o hypercube.o
	ar  rcs libPopGenLib.a  haploid_gt_dis.o haploid_population.o hypercube_function.o hypercube.o


haploid_gt_dis.o:	popgen.h haploid_gt_dis.cpp
	g++ ${CFLAGS} haploid_gt_dis.cpp

haploid_population.o:	popgen.h haploid_population.cpp
	g++ ${CFLAGS} haploid_population.cpp

hypercube_function.o: popgen.h hypercube_function.cpp
	g++ ${CFLAGS} hypercube_function.cpp

hypercube.o: popgen.h hypercube.cpp
	g++ ${CFLAGS} hypercube.cpp
