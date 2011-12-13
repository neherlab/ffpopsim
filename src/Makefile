CXX = g++
CXXFLAGS= -O2 -c -g3 -fPIC


PopGenLib: haploid_gt_dis.o haploid_population.o haploid_clone.o hypercube_function.o hypercube.o sample.o
	ar rcs libPopGenLib.a haploid_gt_dis.o haploid_population.o haploid_clone.o hypercube_function.o hypercube.o sample.o

sample.o: sample.cpp popgen.h
	$(CXX) ${CXXFLAGS} sample.cpp

haploid_gt_dis.o: popgen.h haploid_gt_dis.cpp
	$(CXX) ${CXXFLAGS} haploid_gt_dis.cpp

haploid_population.o: popgen.h haploid_population.cpp
	$(CXX) ${CXXFLAGS} haploid_population.cpp

haploid_clone.o: popgen.h haploid_clone.cpp
	$(CXX) ${CXXFLAGS} haploid_clone.cpp

hypercube_function.o: popgen.h hypercube_function.cpp
	$(CXX) ${CXXFLAGS} hypercube_function.cpp

hypercube.o: popgen.h hypercube.cpp
	$(CXX) ${CXXFLAGS} hypercube.cpp

clean:
	rm -rf *.a *.o
