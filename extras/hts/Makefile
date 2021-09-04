opt=-O3 -g
CXX=icpc
CXXFLAGS=$(opt) -Wall -pedantic -fopenmp -std=c++14 -DNO_BLAS -DUSE_COMPLEX
#                                         ^ for initial build convenience; use BLAS in practice
LDFLAGS=-fopenmp -lufget 
# If using MKL:
# CXXFLAGS=-DUSE_MKL
# LDFLAGS=-mkl
# On platforms supporting avx, include, e.g., -xcore-avx512

SOURCES=hts_test.cpp matrix.cpp mmio.cpp
OBJECTS=$(SOURCES:.cpp=.o)

.cpp.o:
	$(CXX) $(CFLAGS) $(CXXFLAGS) -c $< -o $@

all: hts_test

hts_test: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o hts_test

clean:
	rm -f *.o hts_test

hts_test.o: hts_test_util.hpp hts.hpp hts_decl.hpp hts_impl.hpp hts_impl_def.hpp common.h matrix.h mmio.h dummyFactorize.h
