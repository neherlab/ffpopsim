LIBS += -L/usr/local/lib -lgsl -lgslcblas -lm
HEADERS += \
    ../src/ffpopsim_highd.h \
    ../src/ffpopsim_generic.h \
    ../src/multiLocusGenealogy.h \
    ../src/multi_population.h

SOURCES += \
    ../src/rootedTree.cpp \
    ../src/multiLocusGenealogy.cpp \
    ../src/hypercube_highd.cpp \
    ../src/haploid_highd.cpp \
    ../src/mainsrc.cpp \
    ../src/multi_population.cpp
