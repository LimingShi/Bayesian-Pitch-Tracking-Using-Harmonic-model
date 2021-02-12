# Mandatory
#GTEST = /usr/lib/
#HDF5 = /usr/lib/x86_64-linux-gnu/hdf5/serial


PATH_TO_GTEST=/usr/local
PATH_TO_FFTW3=/usr/local/Cellar/fftw/3.3.8_1
PATH_TO_HDF5=/usr/local/Cellar/hdf5/1.10.5_1
# PATH_TO_eigen=/usr/local/Cellar/eigen/3.3.7







PATH_TO_GTEST_h= -I$(PATH_TO_GTEST)/include/eigen3
PATH_TO_FFTW3_h=-I$(PATH_TO_FFTW3)/include
PATH_TO_HDF5_h=-I$(PATH_TO_HDF5)/include
# PATH_TO_eigen_h= -I$(PATH_TO_GTEST)/include

PATH_TO_GTEST_lib= -L$(PATH_TO_GTEST)/lib
PATH_TO_FFTW3_lib=-L$(PATH_TO_FFTW3)/lib
PATH_TO_HDF5_lib=-L$(PATH_TO_HDF5)/lib





# EXTLIB2=-L$(fftwlib)
# EXTLIB = -L$(HDF5)
# INCLUDE = -I$(HDF5)/include
# INCLUDE2 = -I$(GTESTinclude)/include
# LIBSHDF = $(EXTLIB) -lhdf5_hl -lhdf5 -lz

# Optional - but remove link option -DMKL
# MKLROOT = /opt/intel
#Dynamic linking
# MKL = -m64 -I$(MKLROOT)/mkl/include -I$(MKLROOT)/mkl/include/fftw -L$(MKLROOT)/mkl/lib/intel64 -L$(MKLROOT)/compilers_and_libraries/linux/lib/intel64_lin -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -DALIGNMENT32

FTYPE = -DDOUBLE

#LINK =  -DMKL $(MKL) -DCHEBYSHEV
#LINK = -DMKL $(MKL)
# LINK = $(EXTLIB2) -lfftw3 -DCHEBYSHEV

all:
# 	g++ -g -fno-inline -Wno-write-strings $(FTYPE) $(INCLUDE) $(INCLUDE2)  ../src/single_pitch.cpp ../src/tools.cpp  ../src/th.cpp ../src/vector.cpp fastNLS.cpp $(GTEST)/libgtest.a $(LIBSHDF) $(LINK) -pthread -o fastNLS
	g++ -g -fno-inline -Wno-write-strings $(FTYPE) $(PATH_TO_GTEST_h) $(PATH_TO_FFTW3_h) $(PATH_TO_HDF5_h) ../src/single_pitch.cpp ../src/tools.cpp  ../src/th.cpp ../src/vector.cpp fastNLS.cpp $(PATH_TO_GTEST_lib) -lgtest $(PATH_TO_FFTW3_lib) -lfftw3 -DCHEBYSHEV $(PATH_TO_HDF5_lib) -lhdf5_hl -lhdf5 -lz -pthread -o fastNLS

