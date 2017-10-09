CCXX = g++

CXXFLAGS = -std=c++11 -O3 -lstdc++ -march=native -mtune=native -lm

INCLUDE_FFTW = -lfftw3
INCLUDE_EIGEN = -I/usr/local/include/eigen3
INCLUDE_GSL = -lgsl -lgslcblas
INCLUDE_GNUPLOT_IO = -I /home/artfin/Desktop/repos/gnuplot-iostream  

LINK_BOOST = -lboost_iostreams -lboost_system -lm

########################################################################

POT_SRC = ../POT/src/
POT_OBJ = ../POT/obj/
POT_INCLUDE = -I ../POT/include/

$(POT_OBJ)ar_he_pes.o: $(POT_SRC)ar_he_pes.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(POT_INCLUDE)

########################################################################

one_dimensional_prop: one_dimensional_prop.cpp $(POT_OBJ)ar_he_pes.o
	$(CCXX) $(CXXFLAGS) $(INCLUDE_FFTW) $(INCLUDE_GNUPLOT_IO) $(LINK_BOOST) $(POT_INCLUDE) $^ -o $@ 



