CCXX = g++

CXXFLAGS = -std=c++11 -O3 -lstdc++ -march=native -mtune=native -lm

INCLUDE_FFTW = -lfftw3
INCLUDE_EIGEN = -I/usr/local/include/eigen3
INCLUDE_GSL = -lgsl -lgslcblas

one_dimensional_prop: one_dimensional_prop.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_FFTW) $^ -o $@



