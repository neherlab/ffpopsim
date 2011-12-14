CXX = g++
CXXFLAGS= -O2 -g3 -fPIC


PopGenLib: haploid_gt_dis.o haploid_population.o haploid_clone.o hypercube_function.o hypercube.o sample.o
	ar rcs libPopGenLib.a $^

%.o: %.cpp popgen.h 
	$(CXX) ${CXXFLAGS} -c $<

clean:
	rm -rf *.a *.o
