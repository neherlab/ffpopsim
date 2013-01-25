/**
 * @file mainpage.cpp
 * @brief Main page for the documentation.
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-19
 *
 * @mainpage FFPopSim Documentation
 *
 * @image html quote.png
 *
 * FFPopSim is a C++/Python library used to simulate evolution of genome populations. This is the home page of the C++ documentation.
 *
 * Two kinds of simulations are supported, depending on their genome size:
 * - few loci (up to ~20): the *lowd classes and functions, mainly `haploid_lowd` and `hypercube_lowd`.
 *   The full space of possible genomes is generated and simulations can be performed with any population size.
 *
 * - many loci (> 20, tested up to several tens of thousands): the *highd classes and functions, mainly `haploid_highd` and `hypercube_highd`.
 *   Only the observed genotypes are monitored. If a new genotype is created either by mutation or by recombination, a more resources are required.
 *   Population sizes up to \f$10^6\f$ can be modeled this way. Furthermore, genealogies of single loci can be tracked, and trees can be exported
 *   to standard formats at the end.
 *
 * In addition, FFPopSim includes a specialized subclass of `haploid_highd` used for simulating evolution of the Human Immunodeficiency Virus (HIV),
 * `hivpopulation`. Viruses have a fixed genome size of 10000, and default constructor values are offered to start reasonable simulations quickly.
 * Furthermore, helper I/O routine are included.
 *
 * @section tests-examples Tests and Examples
 *
 * Examples and tests are found in the `tests` folder. They include a few functions that show that the Fast Fourier Transform-based algorithm for
 * recombination yields the same result as the more expensive naive algorithm.
 *
 * @section compiling Compiling and Linking against FFPopSim
 *
 * To build the C++ part of FFPopSim, refer to the INSTALL file.
 * Building your own programs using FFPopSim is then a two-step process, @em compiling @em and @em linking @em.
 * For the compile step, you need to import the header(s) you are using, like in the example above in this page. The FFPopSim headers include recursively all sublibraries such as GSL. Moreover, the compiler must be told where those headers are located; a typical compile call using g++ is the following:
 * \verbatim
 g++ -c -I../pkg/include/ -o example_lowd.o example_lowd.cpp
 \endverbatim
 *
 * For the link step, the following libraries must be used: FFPopSim, gsl, gslcblas. The linker must be told where to find them. Since FFPopSim is compiled statically, its archive must be included in the list of target files. A typical link call using g++ is the following:
 * \verbatim
 g++ example_lowd.o ../pkg/lib/libFFPopSim.a -L../pkg/lib/ -lFFPopSim -lgsl -lgslcblas -o example_lowd
 \endverbatim
 *
 * Of course, compile and link options such as optimization levels can be normally used with FFPopSim.
 *
 * @section python Python Interface
 *
 * FFPopSim also includes bindings to Python 2.7. If you are starting using the library and have Python and all other requirements installed, it is
 * easiest to start out there. Once you are feeling confident with the interface of FFPopSim, the C++ documentation can enlighten you on the internals
 * and details.
 */

