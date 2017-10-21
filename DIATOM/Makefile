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

DIP_SRC = ../DIP/src/
DIP_OBJ = ../DIP/obj/
DIP_INCLUDE = -I ../DIP/include/

$(DIP_OBJ)ar_he_dip.o: $(DIP_SRC)ar_he_dip.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(DIP_INCLUDE)

########################################################################

one_dimensional_prop: one_dimensional_prop.cpp $(POT_OBJ)ar_he_pes.o $(DIP_OBJ)ar_he_dip.o
	$(CCXX) $(CXXFLAGS) $(INCLUDE_FFTW) $(INCLUDE_GNUPLOT_IO) $(LINK_BOOST) $(POT_INCLUDE) $(DIP_INCLUDE) $^ -o $@ 

make_spectrum: make_spectrum.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_GNUPLOT_IO) $(LINK_BOOST) $(INCLUDE_FFTW) $^ -o $@

